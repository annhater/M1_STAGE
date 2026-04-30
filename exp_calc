#IMPORTS
import os
import MDAnalysis as mda
from MDAnalysis.analysis import distances, rms, align
import matplotlib.pyplot as plt
import numpy as np
import warnings
import pandas as pd

# SETUP
# set working dir
rootdir = '/home/sdv/m1isdd/aperova/Documents/M1_STAGE/Data/simulations_1HSI'
os.chdir(rootdir)
warnings.filterwarnings('ignore')

#FUNCTIONS
def dist_50_50(PDB_file):
    u = mda.Universe(PDB_file)
    chA_50 = u.select_atoms('name CA and resid 50')
    chB_50 = u.select_atoms('name CA and resid 149')
    dist = distances.distance_array(chA_50, chB_50,
                                    box=u.dimensions)
    dist = np.float64(dist[0][0])
    sim_name = os.path.basename(os.path.normpath(PDB_file)).replace("_final_renum.pdb", "")
    print(f'Distance 50-50 for {sim_name}: {dist:.3f} Angstrom')
    return dist

def RMSD_flaps(PDB_file, exp_file):
    u = mda.Universe(PDB_file)
    sim_name = os.path.basename(os.path.normpath(PDB_file)).replace("_final_renum.pdb", "")
    ref = mda.Universe(exp_file)
    flap_selection = 'backbone and (resid 43-58 or resid 142-157)'
    #flap_chA = 'backbone and resid 43-58'
    #flap_chB = 'backbone and resid 142-157'
    
    R = rms.RMSD(u,  # universe to align
            ref,  # reference universe or atomgroup
            select=flap_selection  # group to superimpose and calculate RMSD
            #groupselections=[flap_chA, flap_chB]  # groups for RMSD
            )  # frame index of the reference
    R.run()
    df = pd.DataFrame(R.results.rmsd,
                columns=['Frame', 'Time (ns)',
                        'RMSD_flaps'])
    print(f'RMSD_flaps of {sim_name} with {os.path.normpath(exp_file).split('/')[-1]}:')
    print(df)

    return R.results.rmsd

def RMSD_tips(PDB_file, exp_file):
    u = mda.Universe(PDB_file)
    sim_name = os.path.basename(os.path.normpath(PDB_file)).replace("_final_renum.pdb", "")
    ref = mda.Universe(exp_file)
    tips_selection = 'backbone and (resid 48-53 or resid 147-152)'
    tip_chA = 'backbone and resid 48-53'
    tip_chB = 'backbone and resid 147-152'
    
    R = rms.RMSD(u,  # universe to align
            ref,  # reference universe or atomgroup
            select=tips_selection  # group to superimpose and calculate RMSD
            #groupselections=[tip_chA, tip_chB]  # groups for RMSD
            )  # frame index of the reference
    R.run()
    df = pd.DataFrame(R.results.rmsd,
                columns=['Frame', 'Time (ns)',
                        'RMSD_tips'])
    print(f'RMSD_tips of {sim_name} with {os.path.normpath(exp_file).split('/')[-1]}:')
    print(df)

    return R.results.rmsd

#MAIN

# 1. Extract final pbds from all sims

cond1_files = [] # For simulations_1HSI_APO_MP
cond2_files = [] # For simulations_1HSI_APO_DP

for subdir, dirs, files in os.walk(rootdir):
    pdb_files_pathnames = [f for f in files if f.endswith("_final_renum.pdb")]

    for pdb_name in pdb_files_pathnames:
        pdb_path = os.path.join(subdir, pdb_name)
        if "simulations_1HSI_APO_MP" in subdir:
            cond1_files.append(pdb_path)
        elif "simulations_1HSI_APO_DP" in subdir:
            cond2_files.append(pdb_path)

# 2. Extract pdbs from experimental
experimental_files = ["experimental/bent_450.pdb",
                    "experimental/extend_380.pdb",
                    "experimental/Inward_bent_219.pdb"]

""" # Calculate distances
for pdb in cond1_files:
    dist_50_50(pdb)

for pdb in cond2_files:
    dist_50_50(pdb)

for pdb in experimental_files:
    dist_50_50(pdb) """


"""     
    rmsd_tips = rms.rmsd(u.select_atoms('resid 43:58 and resid 142:157').positions,  # coordinates to align
        ref.select_atoms('resid 43:58 and resid 142:157').positions,  # reference coordinates
        center=True,  # subtract the center of geometry
        superposition=True)  # superimpose coordinates """

# 2. Calculate RMSD
for pdb in cond1_files:
    rmsd_flaps_dico = {}
    for exp_file in experimental_files:
        #print(exp_file[i], type(exp_file[i]))
        #RMSD_tips(pdb, exp_file)
        rmsd_flaps_dico[exp_file] = RMSD_flaps(pdb, exp_file)
    #print(rmsd_tips_dico)

for pdb in cond1_files:
    rmsd_tips_dico = {}
    for exp_file in experimental_files:
        #print(exp_file[i], type(exp_file[i]))
        #RMSD_tips(pdb, exp_file)
        rmsd_tips_dico[exp_file] = RMSD_tips(pdb, exp_file)
    #print(rmsd_tips_dico)

for pdb in cond2_files:
    rmsd_flaps_dico = {}
    for exp_file in experimental_files:
        #print(exp_file[i], type(exp_file[i]))
        #RMSD_tips(pdb, exp_file)
        rmsd_flaps_dico[exp_file] = RMSD_flaps(pdb, exp_file)
    #print(rmsd_tips_dico)

for pdb in cond2_files:
    rmsd_tips_dico = {}
    for exp_file in experimental_files:
        #print(exp_file[i], type(exp_file[i]))
        #RMSD_tips(pdb, exp_file)
        rmsd_tips_dico[exp_file] = RMSD_tips(pdb, exp_file)
    #print(rmsd_tips_dico)

""" 
# 2.2. RMSD for c1
for pdb in cond1_files:
    resids, rmsd = avg_rmsd(gro, xtc)
    all_resids_c1.append(resids)
    all_rmsd_c1.append(rmsd)

# 2.3. RMSD for c2
for gro, xtc in cond2_files:
    resids, rmsd = avg_rmsd(gro, xtc)
    all_resids_c2.append(resids)
    all_rmsd_c2.append(rmsd) """
