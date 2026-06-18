
#FUNCTIONS
# a. Open interactions file
def find_hbonds(df):  
    hydrogen_bonds = [col for col in df.columns if 'hb' in col]
    return(hydrogen_bonds)

# b. Analyse hbonds
def hbond_freq(interactions_file):
    df = pd.read_csv(interactions_file)
    n_frames = len(df)
    if simulation_name in selected_frames_dico:
        frame_range = selected_frames_dico[simulation_name]

    sim_data = {}
    for col in hbonds:
        # Calculate frequency: % of frames where the bond exists (val > 0)
        frequency = (df[col] > 0).sum() / n_frames * 100
        sim_data[col] = frequency
    return sim_data


#MAIN

#interaction_files = ['res_V1.csv', 'res_V11.csv', 'res_V12.csv', 'res_V7.csv', 'res_V8.csv', 'res_V21.csv']
interaction_files = ['res_V12.csv', 'res_V7.csv']
#simulation_names = ['V1', 'V11', 'V12', 'V7', 'V8', 'V21']
simulation_names = ['V12', 'V7']
all_results = []

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

phases_list = ['V12_ph1','V7_ph1','V12_ph2','V7_ph2']
sim_phase_data = {}
for i, sim_file in enumerate(interaction_files):
    simulation_name = simulation_names[i]
    if simulation_name in selected_frames_dico:
        frame_range = selected_frames_dico[simulation_name]
        df = pd.read_csv(sim_file)
        df_trans_state = df.loc[frame_range]
        #n_frames_ts = len(df_trans_state)
        
    if simulation_name in phase_limit_dico:
        phase_limit = phase_limit_dico[simulation_name]
    for i,sim_phase in enumerate(phases_list): 
        #add phase one of simulation x to dico
        sim_phase_data[i] = df_trans_state.loc[:phase_limit]
        #add phase two of simulation to dico
        sim_phase_data[i+1] = df_trans_state.loc[phase_limit:]
    sim_data = {}
    for sim_phase in sim_phase_data:
        hbonds = find_hbonds(sim_phase)
    for col in hbonds:
        # Calculate frequency: % of frames where the bond exists (val > 0)
        frequency = (df_trans_state[col] > 0).sum() / n_frames_ts * 100
        sim_data[col] = frequency
    all_results.append(sim_data)
freq_table = pd.DataFrame(all_results)
freq_table.index = phases_list
filtered_cols = freq_table.columns[(freq_table > 0).any()]
filtered_freq_table = freq_table[filtered_cols]
filtered_freq_table = filtered_freq_table.T

filtered_freq_table.to_csv('/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Manips/Tables/interactions_V7_V12.csv', index=True, header=True, decimal=".", float_format="%.5f")
