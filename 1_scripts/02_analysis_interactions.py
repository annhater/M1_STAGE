# Code to search for interactions (hydrogen bonds) present in simulations;
# Specifically around ASP_25 (chain A) or ASP_124 (chain B), plus nearby residues
# Goal: see the differences in interactions of two conditions Monoprotonated (MP) and Unprotoated (DP) 


#IMPORTS
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

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
    resid_list = ['23','24','25','26','27','28','123','124','125','126','128']
    # Filter columns for H-bonds involving Asp 25/124 and near residues
    hb_cols = [col for col in df.columns if 'hb' in col and any(resid in col for resid in resid_list)]
        
    sim_data = {}
    for col in hb_cols:
        # Calculate frequency: % of frames where the bond exists (val > 0)
        frequency = (df[col] > 0).sum() / n_frames * 100
        sim_data[col] = frequency
        
    return sim_data

#MAIN
"""
concat_df = pd.concat(interactions_db)

lengths = [1001, 1001, 1001, 1001, 501, 501]
labels  = ['V1', 'V11', 'V12', 'V7', 'V8', 'V21']

concat_df['simulation'] = ''

start = 0
for L, lab in zip(lengths, labels):
    end = start + L
    concat_df.iloc[start:end, concat_df.columns.get_loc('simulation')] = lab
    start = end

freq_table = pd.DataFrame()

concat_df = concat_df.set_index('simulation', append=False)   # or append=True to keep original integer index
freq_table['hbss_ASP_124_OD2_ASP_25_OD1'] = concat_df['hbss_ASP_124_OD2_ASP_25_OD1'].groupby('simulation').sum()
freq_table['hbss_ASP_124_OD2_ASP_25_OD2'] = concat_df['hbss_ASP_124_OD2_ASP_25_OD2'].groupby('simulation').sum()
freq_table['hbsb_ASP_25_OD1_GLY_126_N'] = concat_df['hbsb_ASP_25_OD1_GLY_126_N'].groupby('simulation').sum()
freq_table['hbsb_ASP_25_OD2_GLY_126_N'] =  concat_df['hbsb_ASP_25_OD2_GLY_126_N'].groupby('simulation').sum()
freq_table['hbsb_ASP_25_OD1_THR_125_N'] = concat_df['hbsb_ASP_25_OD1_THR_125_N'].groupby('simulation').sum()
freq_table['hbsb_ASP_25_OD2_THR_125_N'] = concat_df['hbsb_ASP_25_OD2_THR_125_N'].groupby('simulation').sum()
freq_table['hbsb_THR_125_O_THR_26_OG1'] = concat_df['hbsb_THR_125_O_THR_26_OG1'].groupby('simulation').sum()
freq_table['hbsb_THR_125_OG1_THR_26_N'] = concat_df['hbsb_THR_125_OG1_THR_26_N'].groupby('simulation').sum()
 """


interaction_files = ['res_V1.csv', 'res_V11.csv', 'res_V12.csv', 'res_V7.csv', 'res_V8.csv', 'res_V21.csv']
simulation_names = ['V1', 'V11', 'V12', 'V7', 'V8', 'V21']
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

fig = plt.figure(figsize=(16, 16))
sns.heatmap(freq_percent_filtered, 
            #annot=True,       # Shows the % values in the boxes
            #fmt=".1f",        # 1 decimal point
            cmap="RdYlGn",
            square=True,
            cbar=True,
            #linewidths=0.3)
            )
#ax.set_aspect(freq_percent.shape[1] / freq_percent.shape[0])
plt.title("Hydrogen Bond Frequency", fontsize=16)
plt.ylabel("Simulation")
plt.xlabel("Interaction Type")
plt.xticks(rotation=45, ha='right', fontsize=8)
plt.tight_layout()

plt.show()
#save plot
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/interactions_heatmap.png", bbox_inches='tight', dpi=300)
