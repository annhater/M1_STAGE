#Explain code here

#we need to calculate d_50, and d_50_179
#dico with criterions? or three diff if statements
#that assign conformations
#calculate over the course of the MD

#plot

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
def cal_d_50_179(distance_df):
    d_50_179 = distance_df["50_179"]
    return d_50_179

phase_limit_dico = {
    'V7': 198,
    'V8': 180,
    'V21': 155,
    'V12': 35,
    'V11': 0,
    'V1': 517,
}


#MAIN

# 1. fetch all distance files
#store all file paths
distance_files_paths = []

for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if file.startswith("distance_"):
	        distance_files_paths.append(os.path.join(subdir, file))

# 2. iterate on each distance file
sim_results = []
all_frequencies = []

# 2.1 create pd DF for each distance file
for i, distance_file in enumerate(distance_files_paths):
    simulation_name = os.path.basename(os.path.normpath(distance_file)).replace("distance_", "").replace(".csv", "")
    #transform into pd dataframe
    distance_df = pd.read_csv(distance_file)
    #transform GROMACS distances in nm to Angstrom (1 nm = 10 A)
    distance_df = distance_df * 10
    distance_df.index.names = ['step']

    # Select only phase 2
    limit = phase_limit_dico.get(simulation_name, 250)
    stable_phase_df = distance_df.iloc[limit:].copy()

    # 2.2 compute d_50-50 & d_50-179 for each distance file
    stable_phase_df['d_50-50'] = cal_d_50(stable_phase_df)
    stable_phase_df['d_50-179'] = cal_d_50_179(stable_phase_df)

    stable_phase_df = stable_phase_df.iloc[:, -2:]
    # 3. Define conformation definition conditions
    conditions = [
    (stable_phase_df['d_50-50'] >= 7) & (stable_phase_df['d_50-179'] >= 4.8),
    (stable_phase_df['d_50-50'] < 7) & (stable_phase_df['d_50-179'] >= 4.8),
    (stable_phase_df['d_50-50'] < 7) & (stable_phase_df['d_50-179'] < 4.8)
        
    ]
    values = ['extended', 'bent', 'bent-inward']
    stable_phase_df['conformation'] = np.select(conditions, values, default='extended')
    # write csv with distances and corresponding conformations
    src_path  = r"/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/"
    p = Path(src_path).parent.joinpath(f"Manips/Tables/{simulation_name}_dist_conf.csv")
    stable_phase_df.to_csv(p, index=True, header=True, decimal=".", float_format="%.3f")


    # 4. Calculate frequencies
    counts = stable_phase_df['conformation'].value_counts(normalize=True) * 100
    counts.name = simulation_name
    all_frequencies.append(counts)

# Combine all into summary table
summary_df = pd.concat(all_frequencies, axis=1).fillna(0).T
summary_df.to_csv('/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Tables/conf_freq.csv', index=True, header=True, decimal=".", float_format="%.3f")

#  Plot
# We will store the percentage of time spent in each conformation
data_for_plot = []

conf_labels = ['Extended', 'Bent', 'Bent-inward']
sim_names = ['V7', 'V8', 'V21', 'V1', 'V11', 'V12']

fig, ax = plt.subplots(figsize=(10, 6))

for freq_series in all_frequencies:
    proportions = [
        freq_series.get('extended', 0),
        freq_series.get('bent', 0),
        freq_series.get('bent-inward', 0)
    ]
    data_for_plot.append(proportions)

# Convert to a format easy for plotting
plot_df = pd.DataFrame(data_for_plot, columns=conf_labels, index=sim_names)

# Plotting a horizontal stacked bar chart
plot_df.plot(kind='barh', stacked=True, ax=ax, color=['lightpink', 'mediumturquoise', 'royalblue'], edgecolor='white')

ax.set_xlabel('Proportion of Simulation Time')
ax.set_title('Conformational Frequency per Simulation')
ax.legend(title="Conformation", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
#save plot
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/conf_freq.png", bbox_inches='tight', dpi=300)


# Updated mapping for plotting
conf_map = {'extended': 2, 'bent': 1, 'bent-inward': 0}
conf_colors = {'extended': 'cyan', 'bent': 'pink', 'bent-inward': 'gray'}

# Determine grid size for subplots
num_files = len(distance_files_paths)
# Adjust ncols and nrows to create a somewhat square grid
num_cols = min(num_files, 3)
num_rows = math.ceil(num_files / num_cols)


fig, axes = plt.subplots(num_rows, num_cols, figsize=(num_cols * 8, num_rows * 4), sharex=True)
axes = axes.flatten()

for i, distance_file in enumerate(distance_files_paths):
    ax = axes[i]
    sim_name = os.path.basename(distance_file).replace("distance_", "").replace(".csv", "")
    
    # Load and calculate (using your existing logic)
    df = pd.read_csv(distance_file) * 10
    df['d_50-50'] = cal_d_50(df)
    df['d_50-179'] = cal_d_50_179(df)
    
    # Define Conformations
    conditions = [
        (df['d_50-50'] >= 7) & (df['d_50-179'] >= 4.8),
        (df['d_50-50'] < 7) & (df['d_50-179'] >= 4.8),
        (df['d_50-50'] < 7) & (df['d_50-179'] < 4.8)
    ]
    df['conformation'] = np.select(conditions, ['extended', 'bent', 'bent-inward'], default='extended')
    
    # Map to numbers for plotting
    df['conf_num'] = df['conformation'].map(conf_map)
    
    # Scatter plot: color by conformation type
    for conf_type, color in conf_colors.items():
        mask = df['conformation'] == conf_type
        ax.scatter(df.index[mask], df['conf_num'][mask], color=color, label=conf_type, s=10)

    # Visualizing the Phase Limit
    limit = phase_limit_dico.get(sim_name, 0)
    ax.axvline(x=limit, color='grey', linestyle=':', alpha=0.6)
    trans = ax.get_xaxis_transform()
    plt.text(limit, .15, 'Phase Limit', transform=trans, color='grey', rotation='vertical')

    ax.set_yticks([0, 1, 2])
    ax.set_yticklabels(['Inward', 'Bent', 'Extended'])
    ax.set_title(f"Conformation Transitions for {sim_name}")
    ax.grid(axis='x', alpha=0.3)

plt.tight_layout()
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/conf_timeline.png", bbox_inches='tight', dpi=300)
