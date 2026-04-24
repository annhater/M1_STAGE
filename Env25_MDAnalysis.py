# This code analyses the environment around CA (within 5A) with MDAnslysis,
# and creates tables with frequency of appearance of one atom within 5A throughout the simulation. 


#IMPORTS
import os
from pathlib import Path
import pandas as pd
import MDAnalysis as mda

#SETUP
rootdir = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/simulations_1HSI/'
os.chdir(rootdir)

#FUNCTIONS
def get_dynamic_environment(gro, xtc, cutoff=5.0):
    u = mda.Universe(gro, xtc)
    # Define your targets: ASP 25 on both chains
    # In MDA, we can use segids or chainIDs instead of atom IDs
    CA_asp_chB = u.select_atoms("name CA and resnum 25 and resname ASP and index 1536:3073")
    
    # We want to find atoms around them
    env_atoms = set()

    #print(f"Analyzing {len(u.trajectory)} frames...")
    for ts in u.trajectory[:]:
        # Find atoms within 'cutoff' of the ASP
        around = u.select_atoms(f"around {cutoff} group target", target=CA_asp_chB)
        
        # Add the unique atom names/numbers found in this frame to our set
        for atm in around.atoms:
            env_atoms.add((atm.name, atm.resname, atm.resnum, atm.index))


    # Sort and print the results
    sorted_env = sorted(list(env_atoms), key=lambda x: x[1])
    return sorted_env

def analyze_env_frequency(gro, xtc, cutoff=5.0):
    u = mda.Universe(gro, xtc)
    CA_asp_chB = u.select_atoms("name CA and resnum 25 and resname ASP and index 1536:3073")
    
    # Dictionary to count how many frames each residue is present
    atoms_counts = {}
    n_frames = 0

    for ts in u.trajectory[::1]:
        n_frames += 1
        # Find residues within cutoff (excluding the Asp itself)
        neighbors = u.select_atoms(f"around {cutoff} group CA_asp_chB", CA_asp_chB=CA_asp_chB).residues
        
        for atm in neighbors.atoms:
            atm_id = f"{atm.index}_{atm.name}_{atm.resname}_{atm.resnum}"
            atoms_counts[atm_id] = atoms_counts.get(atm_id, 0) + 1
        

    # Convert to percentage
    freq_data = [{"Residue": k, "Frequency": (v / n_frames) * 100} for k, v in atoms_counts.items()]
    return pd.DataFrame(freq_data).sort_values(by="Frequency", ascending=False)


# MAIN

# 0. Check
env = get_dynamic_environment("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/simulations_1HSI/simulations_1HSI_APO_DP/V1/V1.gro", "/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/simulations_1HSI/simulations_1HSI_APO_DP/V1/md_V1.xtc")
freq_env = analyze_env_frequency("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/simulations_1HSI/simulations_1HSI_APO_DP/V1/V1.gro", "/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/simulations_1HSI/simulations_1HSI_APO_DP/V1/md_V1.xtc")
print("Interesting atoms within 5A:", pd.DataFrame(env))
print("Frequency of atoms appearing within 5A:", pd.DataFrame(freq_env))



# 1. Fetch files from all simulations
sim_files = [] # For all simulations

for subdir, dirs, files in os.walk(rootdir):
    gro_files_pathnames = [f for f in files if f.endswith(".gro")]
    xtc_files_pathnames = [f for f in files if f.endswith(".xtc")]
    for gro_name in gro_files_pathnames:
        sim_name_match = gro_name.replace(".gro", "")
        found_xtc = None
        for xtc_name in xtc_files_pathnames:
            if sim_name_match in xtc_name  and xtc_name.endswith('.xtc'):
                found_xtc = xtc_name
                break
        if found_xtc:
            gro_path = os.path.join(subdir, gro_name)
            xtc_path = os.path.join(subdir, found_xtc)
            sim_files.append((gro_path, xtc_path))

# 1.1. sort file pairs
sim_files.sort()
all_sim_dist = []
for gro, xtc in sim_files:
    sim_name = os.path.basename(os.path.dirname(gro))
    Asp_env_df = pd.DataFrame(get_dynamic_environment(gro, xtc), columns=['atomname', 'resname', 'resid', 'atomic_index'])
    Asp_env_df = Asp_env_df.sort_values(by=['atomic_index']).reset_index(drop=True)

    # 2.7. Write csv with ASP env residues
    src_path = r"/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/"
    p = Path(src_path).parent.joinpath(f"Manips/Tables/{sim_name}_ASP_env_atoms.csv")
    Asp_env_df.to_csv(p, index=True, header=True, decimal=".", float_format="  %.3f")
    # 2.8. Calculate frequencies
    Asp_env_freq_df = analyze_env_frequency(gro, xtc)
    p = Path(src_path).parent.joinpath(f"Manips/Tables/{sim_name}_ASP_env_atoms_frequency.csv")
    Asp_env_freq_df.to_csv(p, index=True, header=True, decimal=".", float_format="  %.3f")


#1. Transform Df to indexes = atomnames, rows = frequencies
#2. Merge all df to that columns = simulations, rows = atomnames, values = frequencies

df = pd.read_csv('V1_ASP_env_atoms_frequency.csv')
df.set_index('Residue', inplace=True)

freq_df_paths = ["V1_ASP_env_atoms_frequency.csv",
                "V11_ASP_env_atoms_frequency.csv",
                "V12_ASP_env_atoms_frequency.csv",
                "V7_ASP_env_atoms_frequency.csv",
                "V8_ASP_env_atoms_frequency.csv",
                "V21_ASP_env_atoms_frequency.csv"]

df_merge = pd.DataFrame()
for sim in freq_df_paths:
    df_sim = pd.read_csv(sim)
    df_sim.set_index('Residue', inplace=True)
    df_merge[sim.replace("_ASP_env_atoms_frequency.csv", "")] = df_sim['Frequency']
    #df_merge = pd.merge(df_list, axis=1)
#print(df_merge)
df_merge.to_csv('merged_CA_ASP_env_atoms_frequency.csv')