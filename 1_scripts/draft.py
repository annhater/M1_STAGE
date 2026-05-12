
#IMPORTS
import os
from networkx import display
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import re
import numpy as np
from sklearn.tree import DecisionTreeClassifier, plot_tree

#SETUP
rootdir = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/interactions/'
os.chdir(rootdir)


#FUNCTIONS
# a. Open interactions file
def find_hbonds(df):
    # Ensure df is a DataFrame before trying to access .columns
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Input to find_hbonds must be a pandas DataFrame")
    hydrogen_bonds = [col for col in df.columns if 'hb' in col]
    return hydrogen_bonds

# b. Analyse hbonds (This function is not used in the refactored MAIN, so it can be removed or kept if needed elsewhere)
def hbond_freq(interactions_file):
    df = pd.read_csv(interactions_file)
    n_frames = len(df)
    # 'simulation_name' and 'hbonds' are not defined in this scope, making this function problematic
    # It seems like this function was intended for overall simulation frequency, not phase-specific
    if 'simulation_name' in globals() and simulation_name in selected_frames_dico:
        frame_range = selected_frames_dico[simulation_name]

    sim_data = {}
    # This part assumes 'hbonds' is globally defined or passed, which is not ideal
    # For this task, phase-specific frequency calculation is handled directly in MAIN
    # for col in hbonds:
    #     frequency = (df[col] > 0).sum() / n_frames * 100
    #     sim_data[col] = frequency
    return sim_data


#MAIN

interaction_files = ['res_V12.csv', 'res_V7.csv']
simulation_names = ['V12', 'V7']

# This dictionary will store all calculated frequencies, keyed by phase name
all_phase_frequencies = {}

selected_frames_dico = {'V7': range(75,250+1),
                    'V12': range(0, 175+1)}

phase_limit_dico = {
    'V7': 198,
    'V8': 180,
    'V21': 155,
    'V12': 35,
    'V11': 1000,
    'V1': 517,
}

# Define the desired order of columns in the final frequency table
phases_list_ordered = ['V12_ph1','V7_ph1','V12_ph2','V7_ph2']

for i, sim_file in enumerate(interaction_files):
    current_simulation_name = simulation_names[i]

    # 1. Load the full DataFrame
    df_full = pd.read_csv(sim_file)

    # 2. Get the transition state slice based on selected_frames_dico
    df_trans_state = pd.DataFrame() # Initialize to avoid UnboundLocalError
    if current_simulation_name in selected_frames_dico:
        frame_range_obj = selected_frames_dico[current_simulation_name]
        # Use df_full.index.isin() for robust selection of indices within the range
        valid_indices = [idx for idx in df_full.index if idx in frame_range_obj]
        df_trans_state = df_full.loc[valid_indices]
    else:
        print(f"Warning: No frame range defined for {current_simulation_name}. Using full DataFrame.")
        df_trans_state = df_full

    if df_trans_state.empty:
        print(f"Warning: Transition state DataFrame is empty for {current_simulation_name}. Skipping.")
        continue

    # 3. Determine phase limit for the current simulation
    phase_limit = phase_limit_dico.get(current_simulation_name, None)

    if phase_limit is None:
        print(f"Warning: No phase limit defined for {current_simulation_name}. Skipping phase separation.")
        continue

    # 4. Separate into phase 1 and phase 2 DataFrames based on the phase_limit
    # Ensure we are slicing based on the *actual indices* present in df_trans_state
    # Get the max index for phase 1 and min index for phase 2 that respect phase_limit
    ph1_indices = df_trans_state.index[df_trans_state.index <= phase_limit]
    ph2_indices = df_trans_state.index[df_trans_state.index > phase_limit]

    df_ph1 = df_trans_state.loc[ph1_indices] if not ph1_indices.empty else pd.DataFrame(columns=df_trans_state.columns)
    df_ph2 = df_trans_state.loc[ph2_indices] if not ph2_indices.empty else pd.DataFrame(columns=df_trans_state.columns)

    n_frames_ph1 = len(df_ph1)
    n_frames_ph2 = len(df_ph2)

    # 5. Find all hydrogen bonds in the entire transition state DataFrame
    # This ensures all phases calculate frequencies for the same set of bonds
    all_hbonds_in_trans_state = find_hbonds(df_trans_state)

    # 6. Calculate frequencies for Phase 1
    ph1_name = f"{current_simulation_name}_ph1"
    current_ph1_frequencies = {}
    for hbond_col in all_hbonds_in_trans_state:
        if hbond_col in df_ph1.columns and n_frames_ph1 > 0:
            frequency = (df_ph1[hbond_col] > 0).sum() / n_frames_ph1 * 100
        else:
            frequency = 0.0 # If hbond_col not in phase 1 or phase 1 is empty
        current_ph1_frequencies[hbond_col] = frequency
    all_phase_frequencies[ph1_name] = current_ph1_frequencies

    # 7. Calculate frequencies for Phase 2
    ph2_name = f"{current_simulation_name}_ph2"
    current_ph2_frequencies = {}
    for hbond_col in all_hbonds_in_trans_state:
        if hbond_col in df_ph2.columns and n_frames_ph2 > 0:
            frequency = (df_ph2[hbond_col] > 0).sum() / n_frames_ph2 * 100
        else:
            frequency = 0.0 # If hbond_col not in phase 2 or phase 2 is empty
        current_ph2_frequencies[hbond_col] = frequency
    all_phase_frequencies[ph2_name] = current_ph2_frequencies

# 8. Convert the dictionary of dictionaries to a DataFrame
# The keys of all_phase_frequencies will become column names, and inner keys will become index
freq_table = pd.DataFrame(all_phase_frequencies)

# 9. Ensure the columns are in the desired order
# Only include columns that were actually generated
existing_ordered_columns = [col for col in phases_list_ordered if col in freq_table.columns]
freq_table = freq_table[existing_ordered_columns]

# 10. Filter out rows (hbonds) where all frequencies across phases are zero
freq_table = freq_table[freq_table.sum(axis=1) > 0]

# Save to CSV
output_path = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Tables/interactions_V7_V12.csv'

# Ensure the directory exists
output_dir = os.path.dirname(output_path)
os.makedirs(output_dir, exist_ok=True)

freq_table.to_csv(output_path, index=True, header=True, decimal=".", float_format="%.5f")

print(f"Frequency table saved to {output_path}")

# Build a decision tree based on the frequencies

#split dataset in features and target variable
feature_cols = freq_table.columns # Assuming the last column is the target variable, adjust if needed
X = freq_table.index # Features
y = freq_table[feature_cols] # Target variable

# Create Decision Tree classifer object
clf = DecisionTreeClassifier()

# Train Decision Tree Classifer
clf = clf.fit(X,y)

#Predict the response for test dataset
y_pred = clf.predict(X_test)
