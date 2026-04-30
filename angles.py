# This code Code for measuring angles throughout the simulation
#1. Import all simulation files
#2. Define function for each angle (psi, phi,)
#3. Assemble all into df with frame/angles, and then into csv file
#4. Create line graph for each simulation for each angle (unite all angles on one graph, indicate angles with colors)

#IMPORTS
import os
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
import mdtraj as md
#import MDAnalysis as mda
from math import pi


#SETUP
rootdir = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/simulations_1HSI/'
os.chdir(rootdir)


#FUNCTIONS
# no functions defined here

#MAIN
#---------------------------------------------------------------------------------------------
# 1. Calculations for phi and psi angles of Asp25 chain A
#---------------------------------------------------------------------------------------------

# 1.0.0. Define residue of interest     
res_indx = 24 #Asp25

# 1.0.1. Setup (empty lists)
all_sim_asp25_phi = []
all_sim_asp25_psi = []
all_angles = []
sim_names = []
all_n_frames = [] # To store n_frame for each trajectory
min_x = []
max_x = []
min_y = []
max_y = []

# 1.0.1. Fetch all simulation files
simulation_files = []

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
            simulation_files.append((pdb_path, xtc_path))

# 1.1. Load trajectories with MDTraj
for pdb, xtc in simulation_files:
    sim_name = os.path.basename(os.path.dirname(pdb))
    sim_names.append(sim_name)
    topology = md.load(pdb).topology
    traj = md.load(xtc,top=topology)
    all_n_frames.append(traj.n_frames) 
    phi_indx,phi_vals = md.compute_phi(traj)
    psi_indx,psi_vals = md.compute_psi(traj)

    asp25_phi = phi_vals[:, res_indx]
    min_asp25_phi, max_asp25_phi = min(asp25_phi),max(asp25_phi)
    min_x.append(min_asp25_phi)
    max_x.append(max_asp25_phi)
    asp25_psi = psi_vals[:, res_indx]
    min_asp25_psi, max_asp25_psi = min(asp25_psi),max(asp25_psi)
    min_y.append(min_asp25_psi)
    max_y.append(max_asp25_psi)
    asp25_angles = [list(x) for x in zip(asp25_phi, asp25_psi)]
    asp25_angles = np.array(asp25_angles)
    all_angles.append(asp25_angles)

# 1.2. Prep for plot
num_simulations = len(sim_names) # Use the actual number of simulations found
num_cols_plot = min(num_simulations, 3)
num_rows_plot = (num_simulations + num_cols_plot - 1) // num_cols_plot

fig, axes = plt.subplots(num_rows_plot, num_cols_plot, figsize=(num_cols_plot * 8, num_rows_plot * 6), squeeze=False)
axes = axes.flatten()
plot_idx = 0
fig.suptitle("Dihedral Map (Phi and Psi) of Asp25 (chain A) for MP (top) and DP (bottom) Simulations", fontsize=22, y=0.95)

lim_x0, lim_x1 = min(min_x), max(max_x)
lim_y0, lim_y1 = min(min_y), max(max_y)


# 1.3. Plot
for i, (current_sim_name, angles_data, n_frames) in enumerate(zip(sim_names, all_angles, all_n_frames)):
    ax = axes[i] # Assign the correct Axes object for plotting
    # Plot scatter points and store the PathCollection object for the colorbar
    sc = ax.scatter(angles_data[:,0], angles_data[:,1], marker="x", c=np.arange(n_frames))
    ax.set_title(f"{current_sim_name}")
    cbar = fig.colorbar(sc, ax=ax) # Create colorbar using the PathCollection and attach to current axes
    cbar.set_label("Frame")
    ax.set_xlabel(r"$\Phi$ Angle [radians]")
    ax.set_ylabel(r"$\Psi$ Angle [radians]")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(lim_x0-0.5, lim_x1+0.5)
    ax.set_ylim(lim_y0-0.5, lim_y1+0.5)
    plot_idx += 1

plt.show()

# 1.4. Save plot
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/Asp_25_phi_psi_angles.png", bbox_inches='tight', dpi=300)

#-------------------------------------------------------------------------------------------
# 2. Calculations for phi and psi angles of Asp25 chain B
#-------------------------------------------------------------------------------------------


# 2.0.0. Define residue of interest     
res_indx = 123 #Asp124

# 2.0.1. Setup (empty lists)
all_sim_asp124_phi = []
all_sim_asp124_psi = []
all_angles = []
#sim_names = []
all_n_frames = [] # To store n_frame for each trajectory
min_x = []
max_x = []
min_y = []
max_y = []

# 2.1. Load trajectories with MDTraj
for pdb, xtc in simulation_files:
    sim_name = os.path.basename(os.path.dirname(pdb))
    topology = md.load(pdb).topology
    traj = md.load(xtc,top=topology)
    all_n_frames.append(traj.n_frames) 
    phi_indx,phi_vals = md.compute_phi(traj)
    psi_indx,psi_vals = md.compute_psi(traj)

    asp124_phi = phi_vals[:, res_indx]
    min_asp124_phi, max_asp124_phi = min(asp124_phi),max(asp124_phi)
    min_x.append(min_asp124_phi)
    max_x.append(max_asp124_phi)
    asp124_psi = psi_vals[:, res_indx]
    min_asp124_psi, max_asp124_psi = min(asp124_psi),max(asp124_psi)
    min_y.append(min_asp124_psi)
    max_y.append(max_asp124_psi)
    asp124_angles = [list(x) for x in zip(asp124_phi, asp124_psi)]
    asp124_angles = np.array(asp124_angles)
    all_angles.append(asp124_angles)

# 2.2. Prep for plot

num_simulations = len(sim_names) # Use the actual number of simulations found
num_cols_plot = min(num_simulations, 3)
num_rows_plot = (num_simulations + num_cols_plot - 1) // num_cols_plot

fig, axes = plt.subplots(num_rows_plot, num_cols_plot, figsize=(num_cols_plot * 8, num_rows_plot * 6), squeeze=False)
axes = axes.flatten()
plot_idx = 0
fig.suptitle("Dihedral Map (Phi and Psi) of Asp25 (chain B) for MP (top) and DP (bottom) Simulations", fontsize=22, y=0.95)

lim_x0, lim_x1 = min(min_x), max(max_x)
lim_y0, lim_y1 = min(min_y), max(max_y)


# 2.3. Plot
for i, (current_sim_name, angles_data, n_frames) in enumerate(zip(sim_names, all_angles, all_n_frames)):
    ax = axes[i] # Assign the correct Axes object for plotting
    # Plot scatter points and store the PathCollection object for the colorbar
    sc = ax.scatter(angles_data[:,0], angles_data[:,1], marker="x", c=np.arange(n_frames))
    ax.set_title(f"{current_sim_name}")
    cbar = fig.colorbar(sc, ax=ax) # Create colorbar using the PathCollection and attach to current axes
    cbar.set_label("Frame")
    ax.set_xlabel(r"$\Phi$ Angle [radians]")
    ax.set_ylabel(r"$\Psi$ Angle [radians]")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(lim_x0-0.5, lim_x1+0.5)
    ax.set_ylim(lim_y0-0.5, lim_y1+0.5)
    plot_idx += 1

plt.show()

# 2.4. Save plot
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/Asp_124_phi_psi_angles.png", bbox_inches='tight', dpi=300)

#-------------------------------------------------------------------------------------------
# 3. Calculations for chi angles of Asp25 chain B
#-------------------------------------------------------------------------------------------


# 3.0.0. Define residue of interest     
res_indx = 123 #Asp124

# 3.0.1. Setup (empty lists)
all_sim_asp124_chi1 = []
all_sim_asp124_chi2 = []
#all_sim_asp124_chi4 = []
#all_sim_asp124_chi3 = []
all_angles = []
#sim_names = []
all_n_frames = [] # To store n_frame for each trajectory
min_x = []
max_x = []
min_y = []
max_y = []

# 3.1. Load trajectories with MDTraj
for pdb, xtc in simulation_files:
    #sim_name = os.path.basename(os.path.dirname(pdb))
    #sim_names.append(sim_name)
    topology = md.load(pdb).topology
    traj = md.load(xtc,top=topology)
    all_n_frames.append(traj.n_frames) 
    chi1_indx,chi1_vals = md.compute_chi1(traj)
    chi2_indx,chi2_vals = md.compute_chi2(traj)

    asp124_chi1 = chi1_vals[:, res_indx]
    min_asp124_chi1, max_asp124_chi1 = min(asp124_chi1),max(asp124_chi1)
    min_x.append(min_asp124_chi1)
    max_x.append(max_asp124_chi1)
    
    asp124_chi2 = chi2_vals[:, res_indx]
    min_asp124_chi2, max_asp124_chi2 = min(asp124_chi2),max(asp124_chi2)
    min_y.append(min_asp124_chi2)
    max_y.append(max_asp124_chi2)

    asp124_angles = [list(x) for x in zip(asp124_chi1, asp124_chi2)]
    asp124_angles = np.array(asp124_angles)
    all_angles.append(asp124_angles)

# 3.2. Prep for plot

num_simulations = len(sim_names) # Use the actual number of simulations found
num_cols_plot = min(num_simulations, 3)
num_rows_plot = (num_simulations + num_cols_plot - 1) // num_cols_plot

fig, axes = plt.subplots(num_rows_plot, num_cols_plot, figsize=(num_cols_plot * 8, num_rows_plot * 6), squeeze=False)
axes = axes.flatten()
plot_idx = 0
fig.suptitle("Dihedral Map (Chi_1 and Chi_2) of Asp25 (chain B) for MP (top) and DP (bottom) Simulations", fontsize=22, y=0.95)

lim_x0, lim_x1 = min(min_x), max(max_x)
lim_y0, lim_y1 = min(min_y), max(max_y)


# 3.3. Plot
for i, (current_sim_name, angles_data, n_frames) in enumerate(zip(sim_names, all_angles, all_n_frames)):
    ax = axes[i] # Assign the correct Axes object for plotting
    # Plot scatter points and store the PathCollection object for the colorbar
    sc = ax.scatter(angles_data[:,0], angles_data[:,1], marker="x", c=np.arange(n_frames))
    ax.set_title(f"{current_sim_name}")
    cbar = fig.colorbar(sc, ax=ax) # Create colorbar using the PathCollection and attach to current axes
    cbar.set_label("Frame")
    ax.set_xlabel(r"$\chi$ Angle [radians]")
    ax.set_ylabel(r"$\chi$ Angle [radians]")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(lim_x0-0.5, lim_x1+0.5)
    ax.set_ylim(lim_y0-0.5, lim_y1+0.5)
    plot_idx += 1

plt.show()

# 3.4. Save plot
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/Asp_124_chi1_chi2_angles.png", bbox_inches='tight', dpi=300)



#-------------------------------------------------------------------------------------------
# 4. Calculations for chi angles of Asp25 chain A
#-------------------------------------------------------------------------------------------


# 4.0.0. Define residue of interest     
res_indx = 24 #Asp25

# 4.0.1. Setup (empty lists)
all_sim_asp25_chi1 = []
all_sim_asp25_chi2 = []
all_angles = []
all_n_frames = [] # To store n_frame for each trajectory
min_x = []
max_x = []
min_y = []
max_y = []

# 4.1. Load trajectories with MDTraj
for pdb, xtc in simulation_files:
    sim_name = os.path.basename(os.path.dirname(pdb))
    topology = md.load(pdb).topology
    traj = md.load(xtc,top=topology)
    all_n_frames.append(traj.n_frames) 
    chi1_indx,chi1_vals = md.compute_chi1(traj)
    chi2_indx,chi2_vals = md.compute_chi2(traj)

    asp25_chi1 = chi1_vals[:, res_indx]
    min_asp25_chi1, max_asp25_chi1 = min(asp25_chi1),max(asp25_chi1)
    min_x.append(min_asp25_chi1)
    max_x.append(max_asp25_chi1)
    
    asp25_chi2 = chi2_vals[:, res_indx]
    min_asp25_chi2, max_asp25_chi2 = min(asp25_chi2),max(asp25_chi2)
    min_y.append(min_asp25_chi2)
    max_y.append(max_asp25_chi2)

    asp25_angles = [list(x) for x in zip(asp25_chi1, asp25_chi2)]
    asp25_angles = np.array(asp25_angles)
    all_angles.append(asp25_angles)

# 4.2. Prep for plot
fig, axes = plt.subplots(num_rows_plot, num_cols_plot, figsize=(num_cols_plot * 8, num_rows_plot * 6), squeeze=False)
axes = axes.flatten()
plot_idx = 0
fig.suptitle("Dihedral Map (Chi_1 and Chi_2) of Asp25 (chain B) for MP (top) and DP (bottom) Simulations", fontsize=22, y=0.95)

lim_x0, lim_x1 = min(min_x), max(max_x)
lim_y0, lim_y1 = min(min_y), max(max_y)


# 4.3. Plot
for i, (current_sim_name, angles_data, n_frames) in enumerate(zip(sim_names, all_angles, all_n_frames)):
    ax = axes[i] # Assign the correct Axes object for plotting
    # Plot scatter points and store the PathCollection object for the colorbar
    sc = ax.scatter(angles_data[:,0], angles_data[:,1], marker="x", c=np.arange(n_frames))
    ax.set_title(f"{current_sim_name}")
    cbar = fig.colorbar(sc, ax=ax) # Create colorbar using the PathCollection and attach to current axes
    cbar.set_label("Frame")
    ax.set_xlabel(r"$\chi$ Angle [radians]")
    ax.set_ylabel(r"$\chi$ Angle [radians]")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(lim_x0-0.5, lim_x1+0.5)
    ax.set_ylim(lim_y0-0.5, lim_y1+0.5)
    plot_idx += 1

plt.show()

# 4.4. Save plot
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/Asp_25_chi1_chi2_angles.png", bbox_inches='tight', dpi=300)