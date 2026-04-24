# Code to calculate RMSD (root mean square deviation) values throughout the MD simulation run
# RMSD is a measurement of the protein structure deviation with respect to the reference structure
# It has been calculated on the backbone atoms, after fitting the trajectory to the starting structure with the backbone atoms.
# Here, we use two different functions: one to calculate RMSD with respect to calculate "average" structure
# to analyse the overall divergence of a frame from average (avg_rmsd);
# And one to see the stepwise evolution of RMSD to see the velocity of conformational movement (stepwise_rmsd).
# The calculated values are used to create plots with matplotlib, and generate statistics csv file


#IMPORTS
import os
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
import matplotlib.pyplot as plt
import numpy as np

# SETUP
# set working dir
rootdir = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/simulations_1HSI/'
os.chdir(rootdir)

# FUNCTIONS
# a. Average RMSD
def avg_rmsd(gro, xtc):
    u = mda.Universe(gro, xtc)
    # 1. alignment of the trajectory
    average = align.AverageStructure(u, u, select='backbone').run()
    ref = average.universe
    aligner = align.AlignTraj(u, ref, select='backbone', in_memory=True).run()

    calphas = u.select_atoms("backbone")
    R = rms.RMSD(calphas, ref.select_atoms("backbone"))
    R.run()
    return calphas.resids, R.results.rmsd

# b. Stepwise RMSD
def stepwise_rmsd(gro, xtc):
    u = mda.Universe(gro, xtc)
    selection = "backbone"
    atoms = u.select_atoms(selection)
    
    stepwise_values = []
    time_points = []
    # 1. Iterate through the trajectory starting from the second frame
    for ts in u.trajectory[1:]:
        # 1.1. Set the 'current' frame as the structure to analyze
        current_coords = atoms.positions.copy()
        # 1.2. Move to the previous frame to set it as the reference
        u.trajectory[ts.frame - 1]
        prev_coords = atoms.positions.copy()
        # 1.3. Move back to current frame to loop 
        u.trajectory[ts.frame]

        # 1.4. Calculate RMSD between t and t-1
        # No alignment is usually done for stepwise unless you want internal change only
        # MDA's rmsd function: rmsd(A, B)
        r = rms.rmsd(current_coords, prev_coords, center=True, superposition=True)
        
        stepwise_values.append(r)
        time_points.append(ts.time)

    # Return in a MDA format: [index, time, rmsd_val]
    # We prepend a 0 for the first frame since there is no t-1 for frame 0
    rmsd_data = np.column_stack([
        np.arange(len(u.trajectory)), 
        u.trajectory.time_points, 
        np.insert(stepwise_values, 0, 0)
    ])
    
    return atoms.resids, rmsd_data

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

# 2. Calculate RMSD

# 2.1. empty lists to store RMSD results
all_rmsd_c1 = []
all_resids_c1 = []
all_rmsd_c2 = []
all_resids_c2 = []

# 2.2. RMSD for c1
for gro, xtc in cond1_files:
    resids, rmsd = avg_rmsd(gro, xtc)
    all_resids_c1.append(resids)
    all_rmsd_c1.append(rmsd)

# 2.3. RMSD for c2
for gro, xtc in cond2_files:
    resids, rmsd = avg_rmsd(gro, xtc)
    all_resids_c2.append(resids)
    all_rmsd_c2.append(rmsd)

# 3. Plot all RMSD (scaled/fixed axes)
# grid size for subplots
num_simulations = len(all_rmsd_c1) + len(all_rmsd_c2)
num_cols_plot = min(num_simulations, 3)
num_rows_plot = (num_simulations + num_cols_plot - 1) // num_cols_plot

fig, axes = plt.subplots(num_rows_plot, num_cols_plot, figsize=(num_cols_plot * 8, num_rows_plot * 6), squeeze=False)
axes = axes.flatten()
fig.suptitle("Average RMSD: Comparison of each frame to the average of the structure", fontsize=22, y=0.92)
plot_idx = 0
line_c1 = None
line_c2 = None

# c1 RMSD plots
for i in range(len(all_rmsd_c1)):
    ax = axes[plot_idx]
    sim_dir_name = os.path.basename(os.path.dirname(cond1_files[i][0]))

    # R.results.rmsd columns: [frame, time, RMSD_total, RMSD_group1, RMSD_group2]
    frames = all_rmsd_c1[i][:, 0] # Extract time
    rmsd_all = all_rmsd_c1[i][:, 2] # RMSD for all

    line_c1, = ax.plot(frames, rmsd_all, color='pink', linewidth=1.5)

    ax.set_ylabel(r"RMSD ($\AA$)")
    ax.set_xlabel("Frame")
    ax.set_title(f"APO_MP_{sim_dir_name}")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.5,7)
    #ax.set_xlim(0, 1000)
    plot_idx += 1

# c2 RMSD plots
for i in range(len(all_rmsd_c2)):
    ax = axes[plot_idx]
    sim_dir_name = os.path.basename(os.path.dirname(cond2_files[i][0]))

    # R.results.rmsd columns: [frame, time, RMSD_total, RMSD_group1, RMSD_group2]
    frames = all_rmsd_c2[i][:, 0] # Extract frame
    rmsd_all = all_rmsd_c2[i][:, 2] # RMSD for all

    line_c2, = ax.plot(frames, rmsd_all, color='cyan', linewidth=1.5)

    ax.set_ylabel(r"RMSD ($\AA$)")
    ax.set_xlabel("Frame")
    ax.set_title(f"APO_DP_{sim_dir_name}")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0.5,7)
    #ax.set_xlim(0.5, 1000)
    plot_idx += 1

# Hide any unused subplots
for j in range(plot_idx, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout(rect=[0, 0, 0.85, 0.92])

# Create the legend

fig.legend(handles=[line_c1, line_c2],     # The line objects
            labels=["Condition 1 (APO_MP)", "Condition 2 (APO_DP)"],   # The labels for each line
            loc="center",  # Position of legend
            bbox_to_anchor=(0.92, 0.5),    # Small spacing around legend box
            title_fontsize='18',
            fontsize='14',
            frameon=True
            )
plt.show()
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/rmsd_average_scaled.png", bbox_inches='tight', dpi=300)

# 4. Plot all RMSD (unscaled, focus on phase 2)
# grid size for subplots
num_simulations = len(all_rmsd_c1) + len(all_rmsd_c2)
num_cols_plot = min(num_simulations, 3)
num_rows_plot = (num_simulations + num_cols_plot - 1) // num_cols_plot

fig, axes = plt.subplots(num_rows_plot, num_cols_plot, figsize=(num_cols_plot * 8, num_rows_plot * 6), squeeze=False)
axes = axes.flatten()
fig.suptitle("Average RMSD: Comparison of each frame to the average of the structure", fontsize=22)
plot_idx = 0
line_c1 = None
line_c2 = None

# c1 RMSD plots
for i in range(len(all_rmsd_c1)):
    ax = axes[plot_idx]
    sim_dir_name = os.path.basename(os.path.dirname(cond1_files[i][0]))

    # R.results.rmsd columns: [frame, time, RMSD_total, RMSD_group1, RMSD_group2]
    frames = all_rmsd_c1[i][:, 0] # Extract time
    rmsd_all = all_rmsd_c1[i][:, 2] # RMSD for all

    line_c1, = ax.plot(frames, rmsd_all, color='pink', linewidth=1.5)

    ax.set_ylabel(r"RMSD ($\AA$)")
    ax.set_xlabel("Frame")
    ax.set_title(f"APO_MP_{sim_dir_name}")
    ax.grid(True, alpha=0.3)
    #ax.set_xlim(0, 1000)
    plot_idx += 1

# c2 RMSD plots
for i in range(len(all_rmsd_c2)):
    ax = axes[plot_idx]
    sim_dir_name = os.path.basename(os.path.dirname(cond2_files[i][0]))

    # R.results.rmsd columns: [frame, time, RMSD_total, RMSD_group1, RMSD_group2]
    frames = all_rmsd_c2[i][:, 0] # Extract frame
    rmsd_all = all_rmsd_c2[i][:, 2] # RMSD for all

    line_c2, = ax.plot(frames, rmsd_all, color='cyan', linewidth=1.5)

    ax.set_ylabel(r"RMSD ($\AA$)")
    ax.set_xlabel("Frame")
    ax.set_title(f"APO_DP_{sim_dir_name}")
    ax.grid(True, alpha=0.3)
    #ax.set_xlim(0.5, 1000)
    plot_idx += 1

# Hide any unused subplots
for j in range(plot_idx, len(axes)):
    fig.delaxes(axes[j])

plt.tight_layout(rect=[0, 0, 0.85, 0.92])

# Create the legend

fig.legend(handles=[line_c1, line_c2],     # The line objects
            labels=["Condition 1 (APO_MP)", "Condition 2 (APO_DP)"],   # The labels for each line
            loc="center",  # Position of legend
            bbox_to_anchor=(0.92, 0.5),    # Small spacing around legend box
            title_fontsize='18',
            fontsize='14',
            frameon=True
            )
plt.show()
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/rmsd_average.png", bbox_inches='tight', dpi=300)


# Below is the code to compare only two simulations (different conditions)

""" 
# load the universe
u = mda.Universe("simulations_1HSI_APO_MP/V7/V7_all_renum.pdb", "simulations_1HSI_APO_MP/V7/md_V7.xtc")
ref = mda.Universe("simulations_1HSI_APO_DP/V1/V1_renum.pdb", "simulations_1HSI_APO_DP/V1/md_V1.xtc")

# FUNCTIONS
# c. RMSD with MDAnalysis
def avg_rmsd(gro, xtc):
    R = rms.RMSD(gro, xtc,
            select="backbone",
            groupselections=["backbone and resid 1-99", "backbone and resid 100-198"])
    R.run()
    return R

# 2. Compute RMSD
R = avg_rmsd(u,ref)

# Plot RMSD
time = R.results.rmsd[:, 1]
rmsd_A = R.results.rmsd[:, 2] # Chain A
rmsd_B = R.results.rmsd[:, 3] # Chain B

fig = plt.figure(figsize=(12, 8))

# Determine grid size for subplots
#fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols * 8, num_rows * 6), squeeze=False)
#axes = axes.flatten()

plt.plot(time, rmsd_A, color = "red", label="Chain A")
plt.plot(time, rmsd_B, color = "blue", label="Chain B")

plt.ylabel(r"RMSD ($\AA$)")
plt.xlabel("Simulation step")
plt.title("RMSD Comparison: Condition 1 (APO_MP_V7) vs Condition 2 (APO_DP_V1) Simulations")
plt.legend(loc='best')
plt.grid(True)
plt.show()
fig.savefig("rmsd.png")
"""