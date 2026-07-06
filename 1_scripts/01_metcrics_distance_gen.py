# Generalised Code for Target Residue Environment Analysis ()
# Change `RES` and `RESNUM` in the Configuration part to run the analysis on any residue.
# UsesMDAnalysis to find all atoms within 5 Å of the target residue's Cα
# Computes per-atom appearance frequency (% of frames) across each simulation
# Saves environment and frequency tables to CSV
# Builds a seaborn heatmap of interaction frequencies

# Needed library: MDAnalysis

# 1. Imports
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda

# 2. Setup & Configuration (residue of interest, paths)

# ── Target residue ───────────────────────────────────────────────────────────
RES    = 'ASP'   # Three-letter residue name
RESNUM = 25      # Residue number in the topology

# ── Atom-index threshold separating chain A from chain B ────────────────────
CHAIN_BOUNDARY = 1536   # atom indices 0..1535 → chain A; 1536+ → chain B
# change according to your system's topology if needed

# ── Paths ────────────────────────────────────────────────────────────────────
ROOTDIR    = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/simulations_1HSI/'
TABLES_DIR = Path('/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Tables/')
FIGURES_DIR= Path('/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/Manips/Figures/')
TABLES_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

SIMULATION_NAMES = ['V1', 'V11', 'V12', 'V7', 'V8', 'V21']
CUTOFF = 5.0

os.chdir(ROOTDIR)


# 3. Helper functions

def get_dynamic_environment(gro: str, xtc: str, resnum: int = RESNUM,
                            cutoff: float = CUTOFF) -> list:
    """
    Return the set of all atoms within `cutoff` Å of the target residue Cα
    across the entire trajectory.
    """
    u = mda.Universe(gro, xtc)
    target_CA = u.select_atoms(f'name CA and resnum {resnum}')
    env_atoms = set()
    for _ in u.trajectory:
        around = u.select_atoms(f'around {cutoff} group target', target=target_CA)
        for atm in around.atoms:
            env_atoms.add((atm.name, atm.resname, atm.resnum, atm.index))
    return sorted(env_atoms, key=lambda x: x[3])


def analyze_env_frequency(gro: str, xtc: str, resnum: int = RESNUM,
                        cutoff: float = CUTOFF) -> pd.DataFrame:
    """
    Compute the % frequency with which each atom appears within `cutoff` Å of
    the target residue Cα.
    """
    u = mda.Universe(gro, xtc)
    target_CA = u.select_atoms(f'name CA and resnum {resnum}')
    atom_counts, n_frames = {}, 0
    for _ in u.trajectory[::1]:
        n_frames += 1
        neighbours = (
            u.select_atoms(f'around {cutoff} group target', target=target_CA).residues
            - target_CA.residues
        )
        for atm in neighbours.atoms:
            key = f'{atm.index}_{atm.name}_{atm.resname}_{atm.resnum}'
            atom_counts[key] = atom_counts.get(key, 0) + 1
    freq_data = [{'Residue': k, 'Frequency': v / n_frames * 100}
                for k, v in atom_counts.items()]
    return pd.DataFrame(freq_data).sort_values('Frequency', ascending=False)


# 4. Discover simulation file pairs
sim_files = []
for subdir, _, files in os.walk(ROOTDIR):
    gro_files = [f for f in files if f.endswith('.gro')]
    xtc_files = [f for f in files if f.endswith('.xtc')]
    for gro_name in gro_files:
        key = gro_name.replace('.gro', '')
        found_xtc = next((x for x in xtc_files if key in x), None)
        if found_xtc:
            sim_files.append((os.path.join(subdir, gro_name),
                                os.path.join(subdir, found_xtc)))
sim_files.sort()
print(f'Found {len(sim_files)} simulation pairs.')

# 5. Compute and save environment/frequency tables
for gro, xtc in sim_files:
    sim_name = Path(gro).parent.name
    print(f'Processing {sim_name} ...')

    env_df = pd.DataFrame(
        get_dynamic_environment(gro, xtc),
        columns=['atomname', 'resname', 'resid', 'atomic_index']
    ).sort_values('atomic_index').reset_index(drop=True)
    env_df.to_csv(TABLES_DIR / f'{sim_name}_{RESNUM}_env_atoms.csv', index=False)

    freq_df = analyze_env_frequency(gro, xtc)
    freq_df.to_csv(TABLES_DIR / f'{sim_name}_{RESNUM}_env_atoms_frequency.csv', index=False)

print('All done.')

# 6. Load frequency tables and build heatmap matrix
all_results = {}
for sim_name in SIMULATION_NAMES:
    freq_path = TABLES_DIR / f'{sim_name}_{RESNUM}_env_atoms_frequency.csv'
    if not freq_path.exists():
        print(f'[WARN] {freq_path} not found.')
        continue

    df_sim = pd.read_csv(freq_path)
    target_label = f'CA_{RES}_{RESNUM}'
    filtered = df_sim[
        df_sim['Residue'].apply(lambda x: 'CA' in x and not x.endswith(target_label))
    ].copy()

    rows = []
    for atom in filtered['Residue'].values:
        parts      = atom.split('_')
        atom_idx   = int(parts[0])
        atom_name  = parts[1]
        res_type   = parts[2]
        res_num    = parts[3]
        chain_lbl  = 'chA' if atom_idx < CHAIN_BOUNDARY else 'chB'
        desc_label = f'{atom_name}_{res_type}{res_num}_{chain_lbl}'
        rows.append({
            'Atom_Identifier': atom,
            'Frequency': float(filtered.loc[filtered['Residue'] == atom, 'Frequency'].iloc[0]),
            'Chain': chain_lbl,
            'Residue_Id': f'{res_type}{res_num}',
            'Descriptive_Label': desc_label
        })
    all_results[sim_name] = pd.DataFrame(rows)

# 7. Assemble heatmap matrix
atom_df      = pd.concat([df.assign(Simulation=n) for n, df in all_results.items()],
                        ignore_index=True)
all_residues = sorted(atom_df['Descriptive_Label'].unique())
heatmap_data = pd.DataFrame(0.0, index=SIMULATION_NAMES, columns=all_residues)

for sim_name, df_sim in all_results.items():
    for label, freq in df_sim.groupby('Descriptive_Label')['Frequency'].max().items():
        heatmap_data.loc[sim_name, label] = freq

heatmap_data = heatmap_data.fillna(0.0)
print(f'Heatmap matrix: {heatmap_data.shape}')

# 8. Plot heatmap
y_labels = [lbl.replace('CA_', '') for lbl in heatmap_data.columns]

fig, ax = plt.subplots(
    figsize=(max(8, len(all_residues) * 0.4), max(4, len(SIMULATION_NAMES) * 1.2))
)
sns.heatmap(
    heatmap_data.T,
    annot=False, cmap='Blues', square=True,
    cbar_kws={'label': 'Frequency (%)'},
    ax=ax
)
ax.set_title(f'Interaction Frequency with {RES}{RESNUM}', fontsize=14)
ax.set_xlabel('Simulation', fontsize=12)
ax.set_ylabel('Residue (within 5 Å)', fontsize=12)
ax.set_yticklabels(y_labels, rotation=0, fontsize=10)
ax.set_xticklabels(SIMULATION_NAMES, rotation=45, ha='right', fontsize=12)
plt.tight_layout()
plt.show()
fig.savefig(FIGURES_DIR / f'env_{RESNUM}_heatmap.png', bbox_inches='tight', dpi=150)
print(f'Saved to {FIGURES_DIR / f"env_{RESNUM}_heatmap.png"}')