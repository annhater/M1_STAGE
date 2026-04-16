# Code to calculate RMSF (root mean square fluctuation) values throughout the MD simulation run
# RMSF is the standard deviation of the atomic positions in the trajectory.
# It has been calculated on the α-carbons, after a fit on the backbone. 
# The calculated values are used to create plots with matplotlib, and generate statistics csv file


#IMPORTS
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align
import matplotlib.pyplot as plt
from MDAnalysis.analysis.rms import RMSF
import os

# SETUP
# set working dir
rootdir = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/simulations_1HSI/'
os.chdir(rootdir)

# FUNCTIONS

# a. RMSF with MDAnalysis
def get_rmsf(gro, xtc):
    u = mda.Universe(gro, xtc)
    # 1. alignment of the trajectory
    average = align.AverageStructure(u, u, select='backbone').run()
    ref = average.universe
    aligner = align.AlignTraj(u, ref, select='backbone', in_memory=True).run()

    # 2. compute RMSF for C-alphas
    calphas = u.select_atoms("name CA")
    rmsfer = RMSF(calphas).run()
    return calphas.resids, rmsfer.results.rmsf

#MAIN

# 1. Fetch files from all simulations
cond1_files = [] # For simulations_1HSI_APO_MP
cond2_files = [] # For simulations_1HSI_APO_DP

for subdir, dirs, files in os.walk(rootdir):
    pdb_files_pathnames = [f for f in files if f.endswith("_renum.pdb")]
    xtc_files_pathnames = [f for f in files if f.endswith(".xtc")]

    for pdb_name in pdb_files_pathnames:
        sim_name_match = pdb_name.replace("_renum.pdb", "").replace("_all", "")

        found_xtc = None
        for xtc_name in xtc_files_pathnames:
            if sim_name_match in xtc_name  and xtc_name.endswith('.xtc'):
                found_xtc = xtc_name
                break
        if found_xtc:
            pdb_path = os.path.join(subdir, pdb_name)
            xtc_path = os.path.join(subdir, found_xtc)

            if "simulations_1HSI_APO_MP" in subdir:
                cond1_files.append((pdb_path, xtc_path))
            elif "simulations_1HSI_APO_DP" in subdir:
                cond2_files.append((pdb_path, xtc_path))

# 1.1. sort file pairs
cond1_files.sort()
cond2_files.sort()

# 2. Calculate RMSF

# 2.1. empty lists to store RMSF results
all_rmsf_c1 = []
all_resids_c1 = []
all_rmsf_c2 = []
all_resids_c2 = []

# 2.2. RMSF for c1
for gro, xtc in cond1_files:
    resids, rmsf = get_rmsf(gro, xtc)
    all_resids_c1.append(resids)
    all_rmsf_c1.append(rmsf)

# 2.3. RMSF for c2
for gro, xtc in cond2_files:
    resids, rmsf = get_rmsf(gro, xtc)
    all_resids_c2.append(resids)
    all_rmsf_c2.append(rmsf)

# 3. Plot all RMSF
fig = plt.figure(figsize=(12, 8))

# 3.1. colors
n_c1 = len(all_rmsf_c1)
n_c2 = len(all_rmsf_c2)
cmap_c1 = plt.get_cmap('PuRd', n_c1 + 1)
cmap_c2 = plt.get_cmap('GnBu', n_c2 + 1)

# c1 RMSF plots
for i in range(len(all_rmsf_c1)):
    sim_dir_name = os.path.basename(os.path.dirname(cond1_files[i][0]))
    plt.plot(all_resids_c1[i], all_rmsf_c1[i], color=cmap_c1(i + 1), alpha=0.7, label=f"APO_MP_{sim_dir_name}")

# c2 RMSF plots
for i in range(len(all_rmsf_c2)):
    sim_dir_name = os.path.basename(os.path.dirname(cond2_files[i][0]))
    plt.plot(all_resids_c2[i], all_rmsf_c2[i], color=cmap_c2(i + 1), alpha=0.7, label=f"APO_DP_{sim_dir_name}")

plt.ylabel(r"RMSF ($\AA$)")
plt.xlabel("Residue")
plt.title("RMSF Comparison: Condition 1 (APO_MP) vs Condition 2 (APO_DP) Simulations")
plt.legend(loc='best')
plt.xlim(-0.5, 200)
plt.ylim(-0, 8.5)
plt.grid(True)
plt.show()
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/rmsf_all_simulations.png")


# Calculate the statistics
summary_stats_dict = {} # dico to store loop results

# Collect simulation names and corresponding RMSF data for individual simulations
all_sim_names = []
all_rmsf_data_combined = []
# For condition 1 individual simulations
for i, (pdb_path, xtc_path) in enumerate(cond1_files):
    # Extract sim_dir_name as used in the plotting cell
    sim_dir_name = os.path.basename(os.path.dirname(pdb_path))
    all_sim_names.append(f"APO_MP_{sim_dir_name}")
    all_rmsf_data_combined.append(all_rmsf_c1[i])

# For condition 2 individual simulations
for i, (pdb_path, xtc_path) in enumerate(cond2_files):
    sim_dir_name = os.path.basename(os.path.dirname(pdb_path))
    all_sim_names.append(f"APO_DP_{sim_dir_name}")
    all_rmsf_data_combined.append(all_rmsf_c2[i])

# Add individual simulation statistics to summary_stats_dict
for sim_name, rmsf_array in zip(all_sim_names, all_rmsf_data_combined):
    summary_stats_dict[sim_name] = pd.Series(rmsf_array).describe()

# Calculate descriptive statistics for all simulations in condition 1 combined
combined_rmsf_c1 = np.concatenate(all_rmsf_c1)
summary_stats_dict['Condition 1'] = pd.Series(combined_rmsf_c1).describe()

# Calculate descriptive statistics for all simulations in condition 2 combined
combined_rmsf_c2 = np.concatenate(all_rmsf_c2)
summary_stats_dict['Condition 2'] = pd.Series(combined_rmsf_c2).describe()

# Create the final DataFrame
summary_stats_df = pd.DataFrame(summary_stats_dict)
#print(summary_stats_df)

# Write csv
summary_stats_df.to_csv('/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/RMSF_stats.csv', index=True, header=True, decimal=".", float_format="%.3f")