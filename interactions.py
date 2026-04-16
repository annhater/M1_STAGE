# Code to search for interactions (hydrogen bonds) present in simulations;
# Specifically around ASP_25 (chain A) or ASP_124 (chain B), plus nearby residues
# Goal: see the differences in interactions of two conditions Monoprotonated (MP) and Unprotoated (DP) 


#IMPORTS
import os
import pandas as pd

#SETUP
rootdir = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/interactions/'
os.chdir(rootdir)

#FUNCTIONS
#add a line to check

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
        
    results = []
    for col in hb_cols:
        # Calculate frequency: % of frames where the bond exists (val > 0)
        frequency = (df[col] > 0).sum() / n_frames * 100
        if frequency > 5.0: # Only care about bonds present > 1% of time
            results.append({'Interaction': col, 'Frequency': frequency})
    
    # Return as a DataFrame sorted by importance
    return pd.DataFrame(results).sort_values(by='Frequency', ascending=False)




#MAIN

# Run for different files to compare

#print("\nV11 (Unstable UP) H-bonds:")
#print(analyze_asp_hbonds('res_V11.csv'))

#print("\nV8 (Stable MP) H-bonds:")
#print(analyze_asp_hbonds('res_V8.csv'))

#print("\nV7 (Stable MP) H-bonds:")
#print(analyze_asp_hbonds('res_V7.csv'))

#print("V12 (Stable UP) H-bonds:")
#print(analyze_asp_hbonds('res_V12.csv'))

#print("\nV21 (Stable MP) H-bonds:")
#print(analyze_asp_hbonds('res_V21.csv'))

#print("\nV1 (Stable MP) H-bonds:")
#analyze_asp_hbonds('res_V1.csv')



# Or make it automatic

# 1. fetch all interaction files
#store all file paths
simulations_list = ['res_V1.csv', 'res_V11.csv', 'res_V12.csv', 'res_V21.csv', 'res_V7.csv', 'res_V8.csv']
interactions_files_paths = []

# 2. iterate on each file
for root, dirs, files in os.walk(rootdir):
    for name in files:
        if name in simulations_list:
            interactions_files_paths.append(os.path.join(rootdir, name))

# 3. use function and print the result
for i, interactions_file in enumerate(interactions_files_paths):
    # get simulation name from filename
    simulation_name = os.path.basename(os.path.normpath(interactions_file)).replace("res_", "").replace(".csv", "")
    print(f"Simulation: {simulation_name}")
    analyze_asp_hbonds(interactions_file)