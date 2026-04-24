# This code explores the calculated distances from "distance_V.csv" files for each simulation
# using defined functions (cal_d_50 and cal_d_50_179) and applies conformation definition filter
# The assigned conformations are used to calculate the frequencies of each conformation during 
# the stable phase of the simulation (phase 2); same logic is applied to visualize the change of conformations
# throughout the MD run on a scatter plot

#IMPORTS
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import math
from pathlib import Path

#SETUP
rootdir = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/simulations_1HSI/'
os.chdir(rootdir)

#FUNCTIONS

# a. d_50 calculation
def cal_d_50 (distance_df):
    d_50 = distance_df["50_149"]
    return d_50

# b. d_50_179 calculation
def cal_d_51_179(distance_df):
    d_51_179 = distance_df["51_179"]
    return d_51_179

phase_limit_dico = {
    'V7': 198,
    'V8': 180,
    'V21': 155,
    'V12': 35,
    'V11': 0, #edited here compared to other scripts to indicate that there was no phase limit (phase 1 never ends) 
    'V1': 517,
}


#MAIN

# 1. Preparation

# 1.1. List to store all file paths
distance_files_paths = []

# 1.2. Fetch all distance files

for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if file.startswith("distance_"):
	        distance_files_paths.append(os.path.join(subdir, file))

# 1.3. Create list to store frequency results
all_frequencies = []

# 2. Iterate on each distance file to calculate frequencies of defined conformations

    # 2.1 Create pd DF for each distance file
for i, distance_file in enumerate(distance_files_paths):
    simulation_name = os.path.basename(os.path.normpath(distance_file)).replace("distance_", "").replace(".csv", "")
    # Transform file into pd dataframe
    distance_df = pd.read_csv(distance_file)
    # Transform GROMACS distances in nm to Angstrom (1 nm = 10 A)
    distance_df = distance_df * 10
    # Add index col name
    distance_df.index.names = ['step']

    # 2.2. Select only phase 2 for frequency calculations
    limit = phase_limit_dico.get(simulation_name, 250) # if not in dico we assume the phase limit at 250
    stable_phase_df = distance_df.iloc[limit:].copy() # copy the df to new df containing only phase 2 steps

    # 2.3. compute d_50-50 & d_50-179 for each distance file
    stable_phase_df['d_50-50'] = cal_d_50(stable_phase_df)
    stable_phase_df['d_51-179'] = cal_d_51_179(stable_phase_df)

    # 2.4. Update df to contain only d_50-50 & d_50-179 columns
    stable_phase_df = stable_phase_df.iloc[:, -2:]
    
    # 2.5. Define different conformations conditions
    conditions = [
    (stable_phase_df['d_50-50'] >= 7) & (stable_phase_df['d_51-179'] >= 4.8),
    (stable_phase_df['d_50-50'] < 7) & (stable_phase_df['d_51-179'] >= 4.8),
    (stable_phase_df['d_50-50'] < 7) & (stable_phase_df['d_51-179'] < 4.8)
        
    ]
    values = ['extended', 'bent', 'bent-inward']

    # 2.6. Add new column to df with corresponding conformations
    stable_phase_df['conformation'] = np.select(conditions, values, default='extended')

    # 2.7. Write csv with distances and corresponding conformations
    src_path  = r"/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/"
    p = Path(src_path).parent.joinpath(f"Manips/Tables/{simulation_name}_dist_conf.csv")
    stable_phase_df.to_csv(p, index=True, header=True, decimal=".", float_format="%.3f")

    # 2.8. Calculate frequencies
    counts = stable_phase_df['conformation'].value_counts(normalize=True) * 100
    counts.name = simulation_name
    all_frequencies.append(counts)

# 3. Make a summary table for frequencies (all simulations: outside of the loop)
summary_df = pd.concat(all_frequencies, axis=1).fillna(0).T
summary_df.to_csv('/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Tables/conf_freq.csv', index=True, header=True, decimal=".", float_format="%.3f")

# 4. Plot
# 4.1. Store the percentage of time spent in each conformation in an empty list
freq_plot = []

# 4.2. List labels for plot
conf_labels = ['Extended', 'Bent', 'Bent-inward']
sim_names = ['V7', 'V8', 'V21', 'V1', 'V11', 'V12']

fig, ax = plt.subplots(figsize=(10, 6))

# 4.3. For each simulation make a list of proportions of the three conformations (three values in the list)
for freq_series in all_frequencies:
    proportions = [
        freq_series.get('extended', 0),
        freq_series.get('bent', 0),
        freq_series.get('bent-inward', 0)
    ]
    freq_plot.append(proportions)

# 4.4. Create a df for plotting
plot_df = pd.DataFrame(freq_plot, columns=conf_labels, index=sim_names)

# 4.5. Plot a horizontal bar chart (stacked)
plot_df.plot(kind='barh', stacked=True, ax=ax, color=['lightpink', 'mediumturquoise', 'royalblue'], edgecolor='white')

# 4.6. Add labels, titles, etc.
ax.set_xlabel('Proportion of Simulation Time')
ax.set_title('Conformational Frequency per Simulation')
ax.legend(title="Conformation", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# 4.7. Save plot
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/conf_freq.png", bbox_inches='tight', dpi=300)

# 5. Visualize the conformations throughout the simulations

# 5.1. Create dictionaries for plotting
conf_choice = {'extended': 2, 'bent': 1, 'bent-inward': 0}
conf_colors = {'extended': 'cyan', 'bent': 'pink', 'bent-inward': 'gray'}

# 5.2. Determine grid size for subplots
num_files = len(distance_files_paths)
num_cols = min(num_files, 3)
num_rows = math.ceil(num_files / num_cols)
fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols * 8, num_rows * 4), sharex=True)
axes = axes.flatten()

# 5.3. Iterate on each file for plotting
for i, distance_file in enumerate(distance_files_paths):
    ax = axes[i]
    sim_name = os.path.basename(distance_file).replace("distance_", "").replace(".csv", "")
    
    # 5.4. Use defined functions for calculations (like previously for freq)
    df = pd.read_csv(distance_file) * 10
    df['d_50-50'] = cal_d_50(df)
    df['d_51-179'] = cal_d_51_179(df)
    
    # 5.3. Define Conformations
    conditions = [
        (df['d_50-50'] >= 7) & (df['d_51-179'] >= 4.8),
        (df['d_50-50'] < 7) & (df['d_51-179'] >= 4.8),
        (df['d_50-50'] < 7) & (df['d_51-179'] < 4.8)
    ]
    values = ['extended', 'bent', 'bent-inward']
    df['conformation'] = np.select(conditions, values, default='extended')
    
    # 5.4. Correspond conformations to mapping choices (0,1,2)
    df['conf_num'] = df['conformation'].map(conf_choice)
    
    # 5.5. Create a scatter plot, color by conformation type
    for conf_type, color in conf_colors.items():
        conf_check = df['conformation'] == conf_type # list booleans to check if the conf in df correspond to dico conf
        # if true: we will plot a point, else we will plot it for different conformation on the next iteration
        ax.scatter(df.index[conf_check], df['conf_num'][conf_check], color=color, label=conf_type, s=10) 
        # only gives the coords(T/F, conf type) for plotting the'true' values, where the checked conformation exists

    # 5.6. Visualize the Phase Limit as a line
    limit = phase_limit_dico.get(sim_name, 0)
    ax.axvline(x=limit, color='grey', linestyle=':', alpha=0.6)
    trans = ax.get_xaxis_transform()
    plt.text(limit, .15, 'Phase Limit', transform=trans, color='grey', rotation='vertical')

    #5.7. Add labels
    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(['Inward', 'Bent', 'Extended'])
    ax.set_title(f"Conformation Transitions for {sim_name}")
    ax.grid(axis='x', alpha=0.3)
plt.tight_layout()

# 5.8. Save plot
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/conf_timeline.png", bbox_inches='tight', dpi=300)
