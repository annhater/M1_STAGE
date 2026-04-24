# Code for calculating d_25 distance throughout the MD simulation run
# we want to extract interactions of Asp25(chain B, i.e. 124) OD2
# with Asp25(chain A) OD1 and OD2 & calculate the distances/change throughout the MD run
# use MDAnalysis for distance calc
# plot with matplotlib, and generate statistics csv file (to do)

#IMPORTS
from matplotlib import colors
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import math
from pathlib import Path
import MDAnalysis as mda
import MDAnalysis.analysis.atomicdistances as ad
from MDAnalysis.analysis import distances
import warnings

#SETUP
rootdir = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/simulations_1HSI/'
os.chdir(rootdir)

# suppress some MDAnalysis warnings when writing PDB files
warnings.filterwarnings('ignore')

#FUNCTIONS

# a. ASP124_OD2_ASP25_OD1_dist
def d_ASP124_OD2_ASP25_OD1(gro, xtc):
    u = mda.Universe(gro, xtc)
    chB_ASP25 = u.select_atoms('name OD2 and resname ASP and resid 25 and index 1538:3074')
    chA_ASP25 = u.select_atoms('name OD1 and resname ASP and resid 25 and index 0:1538')
    
    dists =[]    
    for ts in u.trajectory:
        dist = distances.distance_array(chB_ASP25.positions,
                                        chA_ASP25.positions,
                                        box=u.dimensions)
        dists.append(dist)
    dists = np.ravel(dists)
    return dists

# b. d_ASP124_OD2_ASP25_OD2
def d_ASP124_OD2_ASP25_OD2(gro, xtc):
    u = mda.Universe(gro, xtc)
    chB_ASP25 = u.select_atoms('name OD2 and resname ASP and resid 25 and index 1536:3071')
    chA_ASP25 = u.select_atoms('name OD2 and resname ASP and resid 25 and index 0:1537')
    
    dists =[]    
    for ts in u.trajectory:
        dist = distances.distance_array(chB_ASP25.positions,
                                        chA_ASP25.positions,
                                        box=u.dimensions)
        dists.append(dist)
    dists = np.ravel(dists)
    return dists


# MAIN

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
    col_names = ['124_OD2-25_OD1', '124_OD2-25_OD2']
    u_temp = mda.Universe(gro, xtc)
    num_steps = len(u_temp.trajectory)
    sim_name = os.path.basename(os.path.dirname(gro))


    # 2. Calculate distances
    ASP_distance_df = pd.DataFrame(columns=col_names)
    ASP_distance_df.index.names = ['step']
    ASP_distance_df['124_OD2-25_OD1'] = d_ASP124_OD2_ASP25_OD1(gro, xtc)
    ASP_distance_df['124_OD2-25_OD2'] = d_ASP124_OD2_ASP25_OD2(gro, xtc)
    all_sim_dist.append((sim_name, ASP_distance_df))

    # 2.7. Write csv with ASP distances
    src_path = r"/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/"
    p = Path(src_path).parent.joinpath(f"Manips/Tables/{sim_name}_ASP_distances.csv")
    ASP_distance_df.to_csv(p, index=True, header=True, decimal=".", float_format="%.3f")


# 3. Plot 
#Prep for plot
num_simulations = 6
num_cols_plot = min(num_simulations, 3)
num_rows_plot = (num_simulations + num_cols_plot - 1) // num_cols_plot

fig, axes = plt.subplots(num_rows_plot, num_cols_plot, figsize=(num_cols_plot * 8, num_rows_plot * 6), squeeze=False)
axes = axes.flatten()
colors = ['pink', 'cyan']
plot_idx = 0
fig.suptitle("d$_{25}$ distance between delta oxygens: Condition 1 vs Condition 2 Simulations", fontsize=22, y=0.95)


for (sim_name, sim), ax in zip(all_sim_dist, axes):
    sim.plot(ax=ax, color=['pink', 'cyan'], linewidth=0.5)
    ax.set_title(sim_name)
    ax.plot(sim.index, sim['124_OD2-25_OD1'], color=colors[0], linewidth=0.5)
    ax.plot(sim.index, sim['124_OD2-25_OD2'], color=colors[1], linewidth=0.5)
    ax.set_ylabel(r"distance ($\AA$)")
    ax.set_xlabel("Frame")
    ax.set_title(sim_name)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(1.5,20)
    ax.get_legend().remove()
    #ax.set_xlim(0, 1000)
    plot_idx += 1

# Hide any unused subplots
for j in range(plot_idx, len(axes)):
    fig.delaxes(axes[j])

# Create the legend
handles = [plt.Line2D([0], [0], color=colors[0], label='124_OD2-25_OD1'),
            plt.Line2D([0], [0], color=colors[1], label='124_OD2-25_OD2')]

fig.legend(handles=handles,     # The line objects
            loc="center",  # Position of legend
            bbox_to_anchor=(0.92, 0.5),    # Small spacing around legend box
            title_fontsize='18',
            fontsize='14',
            frameon=True
            )
plt.tight_layout(rect=[0, 0, 0.85, 0.92])
plt.show()

fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/asp_oxygen_distances.png", bbox_inches='tight', dpi=300)

