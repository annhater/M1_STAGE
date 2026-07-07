# -*- coding: utf-8 -*-
"""
Pipeline: Hydrogen Bond Frequency Analysis & Phase Clustering
Simulations: V7, V12 (extensible to more)

Steps:
    1.  Load interaction CSVs
    2.  Slice transition-state frames
    3.  Split into phase 1 / phase 2 by phase_limit
    4.  Compute per-phase H-bond frequencies
    5.  Cluster phases by high-variance bonds (dendrogram, raw freq)
    5b. Cluster phases by Δ-frequency (ph2 − ph1)
    5c. Cluster bonds by frequency profile across phases
    6.  Inter-chain interaction analysis (chain A: res 1–99, chain B: res 100–198)
    7.  PCA across all components – grid of PC pair plots to find phase separation
    8.  Decision Tree for interpretability
"""

# ── Imports ────────────────────────────────────────────────────────────────────
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px

from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier, export_graphviz, plot_tree
from io import StringIO

# ── Configuration ──────────────────────────────────────────────────────────────
ROOT_DIR = '/content/drive/MyDrive/M1_STAGE/Data/interactions/'   # adapt to your env
OUTPUT_DIR = '/content/drive/MyDrive/M1_STAGE/Manips/Tables/'     # adapt to your env
os.makedirs(OUTPUT_DIR, exist_ok=True)

INTERACTION_FILES = ['res_V12.csv', 'res_V7.csv']
SIMULATION_NAMES  = ['V12', 'V7']

# Frames of interest per simulation (transition-state window)
SELECTED_FRAMES = {
    'V12': range(0,   175 + 1),
    'V7':  range(75,  250 + 1),
    # 'V1':  range(412, 587 + 1),
    # 'V8':  range(62,  237 + 1),
    # 'V21': range(62,  237 + 1),
    # 'V11': range(698, 873 + 1),
}

# Index threshold separating phase 1 from phase 2
PHASE_LIMITS = {
    'V12': 35,
    'V7':  198,
    # 'V1':  517,
    # 'V8':  180,
    # 'V21': 155,
    # 'V11': 1000,
}

# Desired column order in the frequency table
PHASES_ORDERED = ['V12_ph1', 'V12_ph2', 'V7_ph1', 'V7_ph2']

# ── Helpers ────────────────────────────────────────────────────────────────────

def find_hbonds(df: pd.DataFrame) -> list[str]:
    """Return column names that correspond to hydrogen bonds (contain 'hb')."""
    return [col for col in df.columns if 'hb' in col]


def phase_frequency(df_phase: pd.DataFrame, hbond_cols: list[str]) -> dict:
    """
    Compute the occurrence frequency (%) of each H-bond in a phase DataFrame.
    Returns 0.0 for empty phases or missing columns.
    """
    n = len(df_phase)
    freqs = {}
    for col in hbond_cols:
        if col in df_phase.columns and n > 0:
            freqs[col] = (df_phase[col] > 0).sum() / n * 100
        else:
            freqs[col] = 0.0
    return freqs


def slice_transition_state(df_full: pd.DataFrame, sim_name: str) -> pd.DataFrame:
    """Return only the transition-state rows defined in SELECTED_FRAMES."""
    if sim_name not in SELECTED_FRAMES:
        print(f"[WARN] No frame range for {sim_name}. Using full DataFrame.")
        return df_full
    valid = [i for i in df_full.index if i in SELECTED_FRAMES[sim_name]]
    return df_full.loc[valid]


def get_residue_numbers(bond_name: str) -> tuple[int, int]:
    """
    Parse a bond column name and return the two residue numbers involved.
    Expected format: ..._RES1_<num1>_..._RES2_<num2>_...
    Falls back to extracting all integers and taking the first two.
    Returns (None, None) if parsing fails.
    """
    import re
    nums = [int(x) for x in re.findall(r'\d+', bond_name)]
    if len(nums) >= 2:
        return nums[0], nums[1]
    return None, None


def classify_chain(resnum: int) -> str | None:
    """Return 'A', 'B', or None based on residue number."""
    if resnum is None:
        return None
    if 1 <= resnum <= 99:
        return 'A'
    if 100 <= resnum <= 198:
        return 'B'
    return None


def is_interchain(bond_name: str) -> bool:
    """Return True if the bond connects one residue in chain A and one in chain B."""
    r1, r2 = get_residue_numbers(bond_name)
    return classify_chain(r1) != classify_chain(r2) and None not in (
        classify_chain(r1), classify_chain(r2))


def high_variance_bonds(freq_table: pd.DataFrame,
                        min_delta: float = 20.0,
                        top_n: int | None = None) -> pd.DataFrame:
    """
    Keep bonds whose frequency changes substantially between at least one ph1/ph2 pair.

    Parameters
    ----------
    freq_table : DataFrame  (bonds × phases)
    min_delta  : minimum |ph2 - ph1| in % for at least one simulation
    top_n      : if set, keep only the top_n bonds by max |delta| across sims
    """
    sim_names = list({c.rsplit('_', 1)[0] for c in freq_table.columns})
    max_delta = pd.Series(0.0, index=freq_table.index)

    for sim in sim_names:
        ph1 = f"{sim}_ph1"
        ph2 = f"{sim}_ph2"
        if ph1 in freq_table.columns and ph2 in freq_table.columns:
            delta = (freq_table[ph2] - freq_table[ph1]).abs()
            max_delta = max_delta.combine(delta, max)

    mask = max_delta >= min_delta
    result = freq_table[mask]
    if top_n is not None:
        result = result.loc[max_delta[mask].nlargest(top_n).index]
    return result


def split_phases(df_ts: pd.DataFrame, sim_name: str):
    """
    Split transition-state DataFrame into phase 1 and phase 2.
    Returns (df_ph1, df_ph2) or raises if no phase limit is defined.
    """
    if sim_name not in PHASE_LIMITS:
        raise ValueError(f"No phase limit defined for {sim_name}.")
    lim = PHASE_LIMITS[sim_name]
    idx = df_ts.index
    ph1 = df_ts.loc[idx[idx <= lim]] if (idx <= lim).any() else pd.DataFrame(columns=df_ts.columns)
    ph2 = df_ts.loc[idx[idx >  lim]] if (idx >  lim).any() else pd.DataFrame(columns=df_ts.columns)
    return ph1, ph2

# ── 1–4: Build frequency table ─────────────────────────────────────────────────
os.chdir(ROOT_DIR)

all_phase_frequencies = {}

for sim_file, sim_name in zip(INTERACTION_FILES, SIMULATION_NAMES):
    df_full   = pd.read_csv(sim_file)
    df_ts     = slice_transition_state(df_full, sim_name)

    if df_ts.empty:
        print(f"[WARN] Empty transition state for {sim_name}. Skipping.")
        continue

    df_ph1, df_ph2 = split_phases(df_ts, sim_name)
    hbonds = find_hbonds(df_ts)          # same bond universe for both phases

    all_phase_frequencies[f"{sim_name}_ph1"] = phase_frequency(df_ph1, hbonds)
    all_phase_frequencies[f"{sim_name}_ph2"] = phase_frequency(df_ph2, hbonds)

    print(f"{sim_name}: {len(df_ph1)} frames in ph1, {len(df_ph2)} frames in ph2, "
          f"{len(hbonds)} H-bonds found.")

# Assemble and order columns
freq_table = pd.DataFrame(all_phase_frequencies)
existing_ordered = [c for c in PHASES_ORDERED if c in freq_table.columns]
freq_table = freq_table[existing_ordered]

# Drop bonds that are zero everywhere
freq_table = freq_table[freq_table.sum(axis=1) > 0]

# Save full table
freq_table.to_csv(os.path.join(OUTPUT_DIR, 'freq_table_all.csv'),
                  index=True, decimal='.', float_format='%.5f')
print(f"\nFrequency table: {freq_table.shape[0]} bonds × {freq_table.shape[1]} phases")

# ── 5: Dendrograms ────────────────────────────────────────────────────────────

# --- 5a: High-variance bonds, cluster PHASES ----------------------------------
# Only keep bonds that change substantially (|Δfreq| ≥ MIN_DELTA_PCT) between
# ph1 and ph2 in at least one simulation.  This removes bonds that are
# consistently present/absent in both phases and therefore carry no transition
# signal, which was causing V12 vs V7 to dominate the clustering.

MIN_DELTA_PCT = 20   # % — tune this; lower = more bonds retained

hv_freq = high_variance_bonds(freq_table, min_delta=MIN_DELTA_PCT)
print(f"Bonds retained after |Δfreq| ≥ {MIN_DELTA_PCT}% filter: {len(hv_freq)}")

if len(hv_freq) < 2:
    print("[WARN] Too few high-variance bonds for clustering. Lower MIN_DELTA_PCT.")
else:
    X_phases = hv_freq.T           # (n_phases, n_bonds)
    y_phases = X_phases.index.tolist()

    linked_phases = linkage(X_phases, method='ward')
    fig, ax = plt.subplots(figsize=(8, 5))
    dendrogram(linked_phases, orientation='top', labels=y_phases,
               distance_sort='descending', ax=ax)
    ax.set_title(f'Phase clustering – high-variance bonds (|Δfreq| ≥ {MIN_DELTA_PCT}%)')
    ax.set_ylabel('Distance')
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'dendrogram_phases_hv.png'), dpi=300)
    plt.show()

# --- 5b: Δ-frequency clustering (ph2 − ph1 per simulation) -------------------
# Removes simulation-identity signal entirely; clusters by the *direction* and
# *magnitude* of change.  If V7 and V12 share a common mechanism they should
# cluster together here even if their absolute frequencies differ.

delta_rows = {}
for sim_name in SIMULATION_NAMES:
    ph1_col, ph2_col = f"{sim_name}_ph1", f"{sim_name}_ph2"
    if ph1_col in freq_table.columns and ph2_col in freq_table.columns:
        delta_rows[sim_name] = freq_table[ph2_col] - freq_table[ph1_col]

delta_df = pd.DataFrame(delta_rows).T
delta_df = delta_df.loc[:, (delta_df != 0).any()]

if len(delta_df) >= 2:
    linked_delta = linkage(delta_df, method='ward')
    fig, ax = plt.subplots(figsize=(8, 5))
    dendrogram(linked_delta, orientation='top', labels=delta_df.index.tolist(),
               distance_sort='descending', ax=ax)
    ax.set_title('Simulation clustering – Δ frequency (ph2 − ph1)\n'
                 'Close = similar transition mechanism')
    ax.set_ylabel('Distance')
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'dendrogram_delta_freq.png'), dpi=300)
    plt.show()

# --- 5c: Cluster BONDS by frequency profile across phases --------------------
# Transpose back so each row is a bond; reveals groups of bonds that change
# together across conditions.

if len(hv_freq) >= 2:
    linked_bonds = linkage(hv_freq, method='ward')
    w = max(10, len(hv_freq) // 3)
    fig, ax = plt.subplots(figsize=(w, 7))
    dendrogram(linked_bonds, orientation='top', labels=hv_freq.index.tolist(),
               distance_sort='descending', leaf_rotation=90, leaf_font_size=7, ax=ax)
    ax.set_title(f'Bond clustering – frequency profile across phases '
                 f'(|Δfreq| ≥ {MIN_DELTA_PCT}%)')
    ax.set_ylabel('Distance')
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'dendrogram_bonds.png'), dpi=300)
    plt.show()

# ── 6: Inter-chain interaction analysis ───────────────────────────────────────
# Chain A: residues 1–99 | Chain B: residues 100–198
# We separate all bonds into intra-A, intra-B, and inter-chain (A↔B),
# then compare their frequencies across phases.

print("\n─── Inter-chain interaction analysis ───")

interchain_mask = pd.Series(
    {bond: is_interchain(bond) for bond in freq_table.index}
)
intraA_mask = pd.Series({
    bond: (classify_chain(get_residue_numbers(bond)[0]) == 'A' and
           classify_chain(get_residue_numbers(bond)[1]) == 'A')
    for bond in freq_table.index
})
intraB_mask = pd.Series({
    bond: (classify_chain(get_residue_numbers(bond)[0]) == 'B' and
           classify_chain(get_residue_numbers(bond)[1]) == 'B')
    for bond in freq_table.index
})

freq_interchain = freq_table[interchain_mask]
freq_intraA     = freq_table[intraA_mask]
freq_intraB     = freq_table[intraB_mask]

print(f"  Intra-chain A bonds : {len(freq_intraA)}")
print(f"  Intra-chain B bonds : {len(freq_intraB)}")
print(f"  Inter-chain A↔B    : {len(freq_interchain)}")

# Save inter-chain table
freq_interchain.to_csv(os.path.join(OUTPUT_DIR, 'freq_interchain.csv'),
                       index=True, decimal='.', float_format='%.5f')

# Heatmap of inter-chain bonds across phases
if not freq_interchain.empty:
    # Apply same high-variance filter so the heatmap stays readable
    ic_hv = high_variance_bonds(freq_interchain, min_delta=MIN_DELTA_PCT)
    if ic_hv.empty:
        ic_hv = freq_interchain   # fall back to all if nothing passes the filter

    fig, ax = plt.subplots(figsize=(max(6, len(ic_hv.columns) * 1.5),
                                    max(4, len(ic_hv) * 0.4)))
    sns.heatmap(ic_hv, annot=True, fmt='.0f', cmap='Blues',
                linewidths=0.3, cbar_kws={'label': 'Frequency (%)'}, ax=ax)
    ax.set_title('Inter-chain (A↔B) H-bond frequencies across phases')
    ax.set_xlabel('Phase')
    ax.set_ylabel('Bond')
    plt.xticks(rotation=30, ha='right')
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'heatmap_interchain.png'), dpi=300,
                bbox_inches='tight')
    plt.show()

    # Δ-frequency bar chart: which inter-chain bonds change most at the transition?
    delta_ic = {}
    for sim_name in SIMULATION_NAMES:
        ph1_col, ph2_col = f"{sim_name}_ph1", f"{sim_name}_ph2"
        if ph1_col in ic_hv.columns and ph2_col in ic_hv.columns:
            delta_ic[sim_name] = ic_hv[ph2_col] - ic_hv[ph1_col]

    if delta_ic:
        delta_ic_df = pd.DataFrame(delta_ic)
        delta_ic_df['mean_delta'] = delta_ic_df.mean(axis=1)
        delta_ic_df = delta_ic_df.sort_values('mean_delta')

        fig, ax = plt.subplots(figsize=(8, max(4, len(delta_ic_df) * 0.4)))
        colors = ['#d73027' if v < 0 else '#1a9850' for v in delta_ic_df['mean_delta']]
        ax.barh(delta_ic_df.index, delta_ic_df['mean_delta'], color=colors)
        ax.axvline(0, color='black', linewidth=0.8)
        ax.set_xlabel('Mean Δ frequency ph2 − ph1 (%)')
        ax.set_title('Inter-chain bond changes at the ph1→ph2 transition\n'
                     'Green = gained, Red = lost')
        plt.tight_layout()
        fig.savefig(os.path.join(OUTPUT_DIR, 'delta_interchain.png'), dpi=300,
                    bbox_inches='tight')
        plt.show()
else:
    print("[INFO] No inter-chain bonds detected. Check residue numbering in column names.")

# ── 7: PCA – all PC pairs ─────────────────────────────────────────────────────
# Fit up to N_PCS components, then display every pair (PC1/2, PC1/3, PC2/3 …)
# as a grid of scatter plots.  Higher PCs sometimes capture phase separation
# better than PC1/PC2 when simulation identity dominates the first axes.

N_PCS = 6   # number of principal components to compute; adjust as needed

phase_csv_map = {
    'V7_ph1':  'V7_ph1_interactions.csv',
    'V7_ph2':  'V7_ph2_interactions.csv',
    'V12_ph1': 'V12_ph1_interactions.csv',
    'V12_ph2': 'V12_ph2_interactions.csv',
}

os.chdir(OUTPUT_DIR)

df_pca_input = pd.DataFrame()
for label, fname in phase_csv_map.items():
    if not os.path.exists(fname):
        print(f"[WARN] {fname} not found, skipping for PCA.")
        continue
    df_tmp  = pd.read_csv(fname)
    hb_cols = find_hbonds(df_tmp)
    df_tmp  = df_tmp[hb_cols].copy()
    df_tmp['label'] = label
    df_pca_input = pd.concat([df_pca_input, df_tmp], ignore_index=True)

if not df_pca_input.empty:
    # Clean: drop zero rows and zero columns
    df_pca_input = df_pca_input.loc[
        (df_pca_input.drop('label', axis=1) != 0).any(axis=1)]
    zero_cols = [c for c in df_pca_input.columns
                 if c != 'label' and df_pca_input[c].sum() == 0]
    df_pca_input.drop(columns=zero_cols, inplace=True)

    X_pca = df_pca_input.drop('label', axis=1)
    y_pca = df_pca_input['label'].reset_index(drop=True)

    n_components = min(N_PCS, X_pca.shape[1], X_pca.shape[0])
    scaler     = StandardScaler()
    X_scaled   = scaler.fit_transform(X_pca)
    pca        = PCA(n_components=n_components)
    components = pca.fit_transform(X_scaled)
    var_expl   = pca.explained_variance_ratio_ * 100

    pc_cols = [f'PC{i+1}' for i in range(n_components)]
    pca_df  = pd.DataFrame(components, columns=pc_cols)
    pca_df['label']      = y_pca.values
    pca_df['simulation'] = pca_df['label'].str.split('_').str[0]
    pca_df['phase']      = pca_df['label'].str.split('_').str[1]

    print(f"\nExplained variance by component:")
    for i, v in enumerate(var_expl):
        print(f"  PC{i+1}: {v:.1f}%  (cumulative: {var_expl[:i+1].sum():.1f}%)")

    # --- Individual interactive scatter: PC1 vs PC2 (quick reference) ---------
    fig = px.scatter(pca_df, x='PC1', y='PC2', color='label',
                     title=f'PCA – PC1 vs PC2  '
                           f'({var_expl[0]:.1f}% + {var_expl[1]:.1f}%)',
                     labels={'PC1': f'PC1 ({var_expl[0]:.1f}%)',
                             'PC2': f'PC2 ({var_expl[1]:.1f}%)'})
    fig.update_traces(marker=dict(size=6, opacity=0.6))
    fig.write_html(os.path.join(OUTPUT_DIR, 'pca_PC1_PC2.html'))
    fig.show()

    # --- Grid of all PC pair plots (static, saved as PNG) ---------------------
    from itertools import combinations
    pairs = list(combinations(pc_cols, 2))
    ncols = 3
    nrows = -(-len(pairs) // ncols)   # ceiling division

    label_colors = {lbl: c for lbl, c in zip(
        pca_df['label'].unique(),
        plt.cm.tab10.colors[:len(pca_df['label'].unique())])}

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 4.5, nrows * 4),
                             squeeze=False)
    for ax_idx, (pcx, pcy) in enumerate(pairs):
        row, col = divmod(ax_idx, ncols)
        ax = axes[row][col]
        for lbl, grp in pca_df.groupby('label'):
            ax.scatter(grp[pcx], grp[pcy],
                       label=lbl, s=12, alpha=0.5,
                       color=label_colors[lbl])
        xi, yi = int(pcx[2:]) - 1, int(pcy[2:]) - 1
        ax.set_xlabel(f'{pcx} ({var_expl[xi]:.1f}%)', fontsize=8)
        ax.set_ylabel(f'{pcy} ({var_expl[yi]:.1f}%)', fontsize=8)
        ax.set_title(f'{pcx} vs {pcy}', fontsize=9)
        ax.tick_params(labelsize=7)

    # Hide unused axes
    for ax_idx in range(len(pairs), nrows * ncols):
        row, col = divmod(ax_idx, ncols)
        axes[row][col].set_visible(False)

    # Single shared legend
    handles, labels_leg = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels_leg, loc='lower right',
               fontsize=8, title='Phase', framealpha=0.8)
    fig.suptitle('PCA grid – all PC pair plots', fontsize=12, y=1.01)
    plt.tight_layout()
    fig.savefig(os.path.join(OUTPUT_DIR, 'pca_grid.png'), dpi=300,
                bbox_inches='tight')
    plt.show()

    # --- Per-PC interactive plots (coloured by phase only, for clarity) --------
    for pcx, pcy in [(f'PC{i}', f'PC{i+1}') for i in range(1, n_components)]:
        xi, yi = int(pcx[2:]) - 1, int(pcy[2:]) - 1
        fig_i = px.scatter(pca_df, x=pcx, y=pcy, color='label',
                           symbol='simulation',
                           title=f'PCA – {pcx} vs {pcy}  '
                                 f'({var_expl[xi]:.1f}% + {var_expl[yi]:.1f}%)',
                           labels={pcx: f'{pcx} ({var_expl[xi]:.1f}%)',
                                   pcy: f'{pcy} ({var_expl[yi]:.1f}%)'})
        fig_i.update_traces(marker=dict(size=6, opacity=0.6))
        fig_i.write_html(os.path.join(OUTPUT_DIR, f'pca_{pcx}_{pcy}.html'))
        fig_i.show()

else:
    print("[INFO] No phase CSVs found; PCA skipped.")

# ── 8: Decision Tree ──────────────────────────────────────────────────────────
# Fit on the aggregated frequency table (phases as samples, bonds as features)
# to find the bonds that best discriminate between sim-phases.

X_dt = X_phases                       # (n_phases, n_bonds)
y_dt = y_phases                       # phase labels

clf = DecisionTreeClassifier(max_depth=4, random_state=42)
clf.fit(X_dt, y_dt)

fig, ax = plt.subplots(figsize=(16, 6))
plot_tree(clf, feature_names=X_phases.columns.tolist(),
          class_names=y_dt, filled=True, rounded=True,
          impurity=False, proportion=False, ax=ax)
ax.set_title('Decision Tree – discriminating bonds between sim-phases')
plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, 'decision_tree.png'), dpi=300, bbox_inches='tight')
plt.show()

print("\nMost discriminating features (top 10 by importance):")
importances = pd.Series(clf.feature_importances_, index=X_phases.columns)
print(importances.sort_values(ascending=False).head(10))
