# M1_STAGE

> All codes created during internship at  IsPP BFA INSERM - CNRS UMR 825, under superision of Dr. Leslie Regad & Marine Baillif at Université Paris Cité from Apr to July of 2026. Topic of the internship: Dynamics of HIV-2 protease (PR2).

![GitHub](https://img.shields.io/badge/GitHub-annhater/M1_STAGE-black?logo=github) ![Build Status](https://img.shields.io/github/actions/workflow/status/annhater/M1_STAGE/ci.yml?branch=main)

## 📋 Table of Contents

- [Features](#features)
- [Repository contents](#repository-contents)

## ℹ️ Project Information

- **👤 Author:** Anna Perova
- **📂 Repository:** [https://github.com/annhater/M1_STAGE](https://github.com/annhater/M1_STAGE)

## General

The present codes allow calculating the following metrics needed to study PR2 interactions:
- For active site stability, check **d<sub>25</sub>** and the codes `01_metrics_d25.py` or `01_metrics_d25.ipynb`.
- For flap dynamic, check **d<sub>50</sub>** and the code `01_metrics_d50.py` or `01_metrics_d50.ipynb`.
- Other inter-residue distances, such as **d<sub>25-50</sub>** for active site - flap distances, can be measured with `01_metrics_distance_gend50_superimposed.R`.
- Evolution of the conformational dynamic over the course of the simulation was measured with **RMSD** using `01_metrics_RMSD.py`.
- Flexible regions of the structure were defined using **RMSF** with `01_metrics_RMSF.py` or `01_metrics_RMSF.ipynb`.
- Changes in the phi and psi angles over the course of the simulation were analysed with `02_analysis_angles.py`.
- Pairwise RMSD and in.py`.
- The projection-based metric **d<sub>50B -> 49A, 50A, 51A</sub>**, used to describe the positioning of flap A over B, is analysed in `02_analysis_teraction` comparisons between simulations were carried out in `02_analysis_pw_RMSD.ipynb` and `02_analysis_pw_interactions.ipynb`.

Additionally, the interactions between the residues of interest were studied with notebooks and scripts dedicated to hydrogen bonds, environment analysis, contact maps, and visualisations.

## Repository contents

### Main project files in `Data`
- `interactions` — contain CSV files that mention all interactions within the simulation.
- `simulations_1HSI` — results of the main simualtions analysed, includes monoprotonated simulations (MP), deprotonated (DP) and simulations with mutations
- `*_contact_matrix_analysis.csv` — contain contacts between residues used to build arc diagrams

### Scripts in `1_scripts`
- `01_metcrics_distance_gen.py` — general code for target-residue environment analysis; it identifies atoms within 5 Å of a selected residue Cα, computes their frequency of appearance across frames, and produces heatmaps and CSV tables.
- `01_metrics_RMSD.py` — code to calculate RMSD throughout the MD run, including average-structure RMSD and stepwise RMSD evolution.
- `01_metrics_RMSD_mut.py` — RMSD analysis workflow adapted for mutant systems.
- `01_metrics_RMSF.py` — code to calculate RMSF on Cα or backbone atoms and export the resulting statistics to CSV.
- `01_metrics_RMSF_mut.py` — RMSF workflow adapted for mutant systems.
- `01_metrics_d25.py` — script to calculate the d25 distance between ASP25 residues on the two chains.
- `01_metrics_d50.py` — script to calculate the d50 distance and classify the flap conformation as extended, bent, or inward-bent.
- `01_metrics_dist_50.R` — R script used for the projection-based flap-positioning analysis.
- `02_analysis_ASP_OD.py` — script to calculate distances involving the ASP25 OD atoms and follow their evolution in time.
- `02_analysis_angles.py` — code to measure phi and psi angles throughout the simulation and export the results.
- `02_analysis_conf_freq.py` — analysis of conformation frequencies based on the calculated distances and stable-phase classification.
- `02_analysis_d50_superimposed.R` — R script for the superimposed flap-positioning analysis.
- `02_analysis_env25.py` — script that analyses the environment around ASP25 within 5 Å and builds frequency tables.
- `02_analysis_experimental.py` — script comparing experimental structures to the MD simulations through RMSD-based metrics.
- `02_analysis_interactions.py` — code to search for hydrogen-bond interactions near ASP25 and ASP124.
- `02_analysis_transition.py` — analysis of interaction frequencies during transition states and clustering of those interactions.
- `03_plot_cont_map.py` — script producing a contact-map plot of residue contacts during the simulation.
- `03_plot_dist.R` — R script used to visualise the distance-based projection analysis.
- `03_plot_interaction_heatmap.py` — code to plot a heatmap of interaction frequencies around ASP25.
- `fix_chains.py` — utility script used to fix chain identifiers in PDB files when needed.

### Notebooks in `2_notebooks`
- `01_metrics_RMSF.ipynb` — notebook version of the RMSF analysis for the different MD conditions.
- `01_metrics_d25.ipynb` — notebook for d25-distance analysis, plotting the trace over time and exporting summary statistics.
- `01_metrics_d50_50_new.ipynb` — notebook for a distance-based analysis of flap-region geometry and phase-dependent behaviour.
- `01_metrics_distance.ipynb` — notebook dedicated to general distance-metric analysis.
- `02_analysis_clustering.ipynb` — notebook for clustering H-bond frequencies in the transition-state window.
- `02_analysis_env25_heatmap.ipynb` — notebook producing a heatmap of the residue environment around ASP25.
- `02_analysis_interactions_clustermap.ipynb` — notebook generating a clustermap of ASP25/124 hydrogen-bond frequencies.
- `02_analysis_interactions.ipynb` — notebook dedicated to the broader interaction analysis workflow.
- `02_analysis_pw_RMSD.ipynb` — notebook for pairwise RMSD analysis between the compared simulations.
- `02_analysis_pw_interactions.ipynb` — notebook for pairwise Jaccard-based interaction analysis between the compared simulations.
- `02_analysis_V12_env25_HB.ipynb` — notebook focusing on ASP25 hydrogen-bond interactions for the V12 system.
- `03_plot_HB_network.ipynb` — notebook for visualising hydrogen-bond networks.
- `03_plot_acr_diagram.ipynb` — notebook producing an arc diagram of intra-atomic interactions.
- `03_plot_contact_map.ipynb` — notebook for plotting residue contact maps.
- `03_plot_env25_chA_graph.ipynb` — notebook building an iGraph network for the ASP25 environment on chain A.
- `03_plot_env25_chB_graph.ipynb` — notebook building an iGraph network for the ASP25 environment on chain B.
- `03_plot_env_heatmap.ipynb` — notebook creating a heatmap of the environment around ASP25.
- `03_plot_ramachandran.ipynb` — notebook producing a Ramachandran plot for ASP25 dihedral angles.


### Report files in `3_report`
- `5_report` contains the files used to produce the report: the LaTeX source, bibliography, and the figures present in the report, together with the raw text and supporting assets.

### Archived files in `4_archive`
- `6_archive` gathers older scripts and notebooks that were not kept in the final workflow or that were later transformed into newer versions of the analysis code.
