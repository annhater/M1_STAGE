
#IMPORTS
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


#SETUP
rootdir = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Tables/'
os.chdir(rootdir)

#MAIN
# 1. Extract frequency data

# 1. Read data

# 1.0. Setup
df_list = []
CA_atoms = []
simulation_names = ['V1', 'V11', 'V12', 'V7']
all_results = {}
residue_dico = {}
# Initialize a set to collect all unique atom identifiers across all simulations
all_unique_atom_identifiers = set()

# 1.1. Open frequency tables, read as df
freq_df_paths = [
    "V1_ASP_chA_env_atoms_frequency.csv",
    "V11_ASP_chA_env_atoms_frequency.csv",
    "V12_ASP_chA_env_atoms_frequency.csv",
    "V7_ASP_chA_env_atoms_frequency.csv"
]
for sim in freq_df_paths:
    sim_name = sim.replace("_ASP_chA_env_atoms_frequency.csv", "")
    df_sim = pd.read_csv(sim)

    # 1.2. Filter only alpha carbons (and that are not Asp25 alpha carbon)
    filtered_CA = df_sim[
        df_sim["Residue"].apply(lambda x: 'CA' in x and '2274_CA_ILE_50' not in x and '2273_CA_ILE_50' not in x)
        ].copy()
#384_CA_LEU_24
#403_CA_ASP_25
#415_CA_THR_26
#429_CA_GLY_27
#436_CA_ALA_28
#1919_CA_LEU_24
#1938_CA_ASP_25
#1951_CA_THR_26
#1965_CA_GLY_27

    # 1.3. Save to new df for graph parameters
    df_new = pd.DataFrame()
    df_new["Atom_Identifier"] = filtered_CA["Residue"] # Keep the full atom identifier
    df_new["Frequency"] = filtered_CA['Frequency']

    # 1.3.1. For each atom, get its residue_id, chain, and atom_name
    chain_label_list = []
    res_id_list = []
    atom_name_list = []


    for atom in filtered_CA["Residue"].values:
        res_str = atom.split('_')
        atom_id = int(res_str[0])
        res_id = '_'.join(res_str[2:])
        atom_name = res_str[1]
        if atom_id < 1536:
            chain_label_list.append("Chain A")
        else:
            chain_label_list.append("Chain B")
        res_id_list.append(res_id)
        # Collect all unique atom identifiers in the global set
        atom_name_list.append(atom_name)
        all_unique_atom_identifiers.add(atom) # Add to the global set

    df_new["Chain"] = chain_label_list
    df_new["Residue_Id"] = res_id_list
    df_new["Atom"] = atom_name_list


    all_results[sim_name] = df_new



# Prepare a list to store parsed atom data
all_atom_data = []

# Iterate through all_results to gather and parse atom identifiers
for sim_name, df_sim_data in all_results.items():
    for _, row in df_sim_data.iterrows():
        atom_identifier = row['Atom_Identifier']
        parts = atom_identifier.split('_')
        atom_id = int(parts[0])
        atom_name = parts[1]
        residue_type = parts[2]
        residue_num = parts[3]

        chain_label = 'chA' if atom_id < 1536 else 'chB'

        # Create a 'descriptive_label' that excludes the atom_id
        descriptive_label = f"{atom_name}_{residue_type}{residue_num}_{chain_label}"

        all_atom_data.append({
            'Atom_Identifier_Full': atom_identifier,
            'Atom_ID': atom_id,
            'Descriptive_Label': descriptive_label,
            'Simulation': sim_name
        })

# Create a DataFrame from the collected data
atom_data_df = pd.DataFrame(all_atom_data)

# Group by Descriptive_Label and check for multiple unique Atom_IDs
duplicate_label_groups = atom_data_df.groupby('Descriptive_Label').filter(lambda x: x['Atom_ID'].nunique() > 1)

# Sort the results for better readability
duplicate_label_groups = duplicate_label_groups.sort_values(by=['Descriptive_Label', 'Atom_ID'])

#print("Atoms with the same descriptive label but different Atom_IDs (potential shifts):")
#display(duplicate_label_groups[['Descriptive_Label', 'Atom_ID', 'Atom_Identifier_Full', 'Simulation']])


# The 'atom_data_df' created in cell 'e0f79337' contains the 'Descriptive_Label' column
# Define the list of unique descriptive labels for the heatmap columns
all_residues = atom_data_df['Descriptive_Label'].unique().tolist()

# Initialize the heatmap_data DataFrame with simulation names as index and merged descriptive labels as columns
heatmap_data = pd.DataFrame(index=simulation_names, columns=all_residues)

# Populate the heatmap_data by aggregating frequencies for each simulation
for sim_name in simulation_names:
    # Get the raw frequency data for the current simulation
    df_sim_raw = all_results[sim_name].copy()

    # Add the Descriptive_Label column to this simulation's DataFrame if it doesn't exist
    # This ensures consistency with how 'all_residues' was generated
    descriptive_labels_for_sim = []
    for atom_identifier in df_sim_raw['Atom_Identifier']:
        parts = atom_identifier.split('_')
        atom_id_val = int(parts[0])
        atom_name = parts[1]
        residue_type = parts[2]
        residue_num = parts[3]
        chain_label = 'chA' if atom_id_val < 1536 else 'chB'
        descriptive_labels_for_sim.append(f"{atom_name}_{residue_type}{residue_num}_{chain_label}")
    df_sim_raw['Descriptive_Label'] = descriptive_labels_for_sim

    # Group by the Descriptive_Label and take the maximum frequency
    merged_frequencies = df_sim_raw.groupby('Descriptive_Label')['Frequency'].max()

    # Assign these merged frequencies to the correct row in heatmap_data
    for label, freq in merged_frequencies.items():
        if label in heatmap_data.columns: # Ensure the label exists in our defined columns
            heatmap_data.loc[sim_name, label] = freq

# Convert to float to handle NaN values properly
heatmap_data = heatmap_data.astype(float)
heatmap_data = heatmap_data.fillna(0)

# Sort the columns (residues) alphabetically for consistent display
#heatmap_data = heatmap_data.reindex(columns=sorted(heatmap_data.columns))


fig, ax = plt.subplots(figsize=(20, 15)) # Significantly increased figure height for larger squares and better visibility

sns.heatmap(heatmap_data.T, # Transpose the DataFrame to swap axes
            annot=False,        # Keep as False for now to avoid clutter with many small cells
            fmt=".1f",          # Format for annotations
            cmap="Blues",      # Color map (green for high frequency, red for low)
            square=True,       # Set to True if you want squares to be equal size
            cbar=True,  # Color of the lines between cells
            ax=ax,              # Pass the ax object to sns.heatmap
            cbar_kws={'label': 'Frequency (%)'}
            )

plt.title("Frequency of Interactions with Ile50 (chain A)", fontsize=16, y=1.05, x=0.66)
plt.xlabel("Simulation", fontsize=12) # Swapped x-axis label
plt.ylabel("Residue (within 5A)", fontsize=12) # Swapped y-axis label

# Use the 'all_residues' (which now contains the unique descriptive labels) directly for y-axis labels
formatted_y_labels = heatmap_data.columns.tolist() # Get the sorted columns from heatmap_data

# Remove 'CA_' prefix from the labels
formatted_y_labels = [label.replace('CA_', '') for label in formatted_y_labels]

plt.yticks(ticks=[x + 0.5 for x in range(len(formatted_y_labels))], labels=formatted_y_labels, rotation=0, fontsize=12) # Apply to y-axis with appropriate rotation, slightly increased fontsize
plt.xticks(ticks=[x + 0.5 for x in range(len(simulation_names))], labels=simulation_names, rotation=45, fontsize=15) # Apply to x-axis, using simulation_names
#plt.tight_layout() # Adjust layout to prevent labels from overlapping
plt.show()
fig.savefig("/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Figures/env_50_chA_heatmap.png", bbox_inches='tight', dpi=300)
