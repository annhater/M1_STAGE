# Code to search for interactions (hydrogen bonds) present in simulations;
# Specifically around ASP_25 (chain A) or ASP_124 (chain B), plus nearby residues
# Goal: see the differences in interactions of two conditions Monoprotonated (MP) and Unprotoated (DP) 

# Created by: Anna Perova

#IMPORTS
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

#SETUP
rootdir = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/interactions/'
os.chdir(rootdir)

#FUNCTIONS
# a. Open interactions file
def find_asp_int(interactions_file):
    interactions_df = pd.read_csv(interactions_file)
    HB_ASP_25_cols = [col for col in interactions_df.columns if 'ASP_25' in col and 'hb' in col]
    for present_interactions in HB_ASP_25_cols:
        if interactions_df[present_interactions].sum() > 0:
            print(present_interactions)

# b. Analyse Asp_bonds
def analyze_asp_hbonds(interactions_file):
    df = pd.read_csv(interactions_file)
    n_frames = len(df)
    resid_list = ['23','24','25','26','27','122','123','124','125','126']
    # Filter columns for H-bonds involving Asp 25/124 and near residues
    hb_cols = [col for col in df.columns if 'hb' in col and any(resid in col for resid in resid_list)]
        
    sim_data = {}
    for col in hb_cols:
        # Calculate frequency: % of frames where the bond exists (val > 0)
        frequency = (df[col] > 0).sum() / n_frames * 100
        sim_data[col] = frequency
        
    return sim_data

# Function to extract the first residue number from your interaction string
def get_res_num(interaction_string):
    # Extracts the first number found in the string
    return int(re.search(r'\d+', interaction_string).group())

#MAIN

interaction_files = ['res_V1.csv', 'res_V11.csv', 'res_V12.csv', 'res_V7.csv']
simulation_names = ['V1', 'V11', 'V12', 'V7']
all_results = []

for sim_file in interaction_files:
    data = analyze_asp_hbonds(sim_file)
    all_results.append(data)
freq_table = pd.DataFrame(all_results)

freq_table.index = simulation_names
freq_table

#lengths_dico = {'V1': 1001, 'V11': 1001, 'V12': 1001, 'V7': 1001, 'V8': 501, 'V21': 501}
#freq_percent = freq_table.div(freq_table.index.map(lengths_dico), axis=0) * 100

# Only keep interactions that occur in at least 25% of frames in at least one simulation
filtered_cols = freq_table.columns[(freq_table > 25).any()]
freq_percent_filtered = freq_table[filtered_cols]

#lengths_dico = {'V1': 1001, 'V11': 1001, 'V12': 1001, 'V7': 1001, 'V8': 501, 'V21': 501}
#freq_percent = freq_table.div(freq_table.index.map(lengths_dico), axis=0) * 100

descriptive_labels = []
for bond in freq_percent_filtered.columns:
    parts = bond.split('_')
    res1 = f"{parts[1]:<3}{parts[2]:>3}_{parts[3]:<3}" 
    res2 = f"{parts[4]:<3}{parts[5]:>3}_{parts[6]:<3}"
    descriptive_labels.append(f"{res2} - {res1}")

for col_name, label in zip(freq_percent_filtered.columns,descriptive_labels):
    freq_percent_filtered = freq_percent_filtered.rename(columns={col_name:label})

sorted_cols = sorted(freq_percent_filtered.columns, key=get_res_num)
df_sorted = freq_percent_filtered[sorted_cols]

#freq_percent_filtered = sorted(freq_percent_filtered.T)
freq_percent_filtered = freq_percent_filtered.reindex(columns=sorted(freq_percent_filtered.columns))
fig = plt.figure(figsize=(16,16))
sns.heatmap(df_sorted.T, 
            #annot=True,       # Shows the % values in the boxes
            #fmt=".1f",        # 1 decimal point
            cmap="Blues",
            square=True,
            cbar=True,
            #linewidths=0.3)
            )
#ax.set_aspect(freq_percent.shape[1] / freq_percent.shape[0])
#plt.title("Fréquence des liaisons hydrogènes", fontsize=16)
plt.xlabel("Simulation")
plt.ylabel("Type d'interaction")
plt.yticks(family='monospace', fontsize=10)
plt.xticks(family='monospace', rotation=90, fontsize=8)
plt.tight_layout()

plt.show()
#save plot
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/interactions_heatmap_2.png", bbox_inches='tight', dpi=150)
