# This is code to analyse interactions and their frequencies that 
# occur during transition states of different MD simulation

# Goal: see the differences in interactions of two conditions Monoprotonated (MP) and Unprotoated (DP)
# Make clustering (HClust in R) to define the rules of interactions

# Created by: Anna Perova


# 1. Test with V12 and V7

import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re

#SETUP
rootdir = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/interactions/'
os.chdir(rootdir)
output_path = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Tables/'

#FUNCTIONS
# a. Open interactions file
def find_hbonds(df):
    # Ensure df is a DataFrame before trying to access .columns
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Input to find_hbonds must be a pandas DataFrame")
    hydrogen_bonds = [col for col in df.columns if 'hb' in col]
    return hydrogen_bonds


#MAIN

interaction_files = ['res_V12.csv', 'res_V7.csv']
simulation_names = ['V12', 'V7']

selected_frames_dico = {'V7': range(162,234+1),
                    'V12': range(0, 72+1)}

phase_limit_dico = {
    'V7': 198,
    'V8': 180,
    'V21': 155,
    'V12': 35,
    'V11': 1000,
    'V1': 517,
}

for i, sim_file in enumerate(interaction_files):
    current_simulation_name = simulation_names[i]

    # 1. Load the full DataFrame
    df_full = pd.read_csv(sim_file)

    # 2. Get the transition state slice based on selected_frames_dico
    df_trans_state = pd.DataFrame()
    if current_simulation_name in selected_frames_dico:
        frame_range_obj = selected_frames_dico[current_simulation_name]
        valid_indices = [idx for idx in df_full.index if idx in frame_range_obj]
        df_trans_state = df_full.loc[valid_indices]

    # 3. Determine phase limit for the current simulation
    phase_limit = phase_limit_dico.get(current_simulation_name)

    # 4. Separate into phase 1 and phase 2 DataFrames based on the phase_limit
    ph1_indices = df_trans_state.index[df_trans_state.index <= phase_limit]
    ph2_indices = df_trans_state.index[df_trans_state.index > phase_limit]

    df_ph1 = df_trans_state.loc[ph1_indices] if not ph1_indices.empty else pd.DataFrame(columns=df_trans_state.columns)
    df_ph2 = df_trans_state.loc[ph2_indices] if not ph2_indices.empty else pd.DataFrame(columns=df_trans_state.columns)

    # 5. Find all hydrogen bonds in the entire transition state DataFrame
    hbonds = find_hbonds(df_trans_state)

    # 6. Interactions for Phase 1
    ph1_name = f"{current_simulation_name}_ph1"
    ph1_presence_absence_df = df_ph1[hbonds].map(lambda x: 1 if x > 0 else 0)
    ph1_output_path = os.path.join(output_path, f"{ph1_name}_interactions.csv")
    ph1_presence_absence_df.to_csv(ph1_output_path, index=True, header=True)

    # 7.Interaction for Phase 2
    ph2_name = f"{current_simulation_name}_ph2"
    ph2_presence_absence_df = df_ph2[hbonds].map(lambda x: 1 if x > 0 else 0)
    ph2_output_path = os.path.join(output_path, f"{ph2_name}_interactions.csv")
    ph2_presence_absence_df.to_csv(ph2_output_path, index=True, header=True)

