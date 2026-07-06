# M1_STAGE

> All codes created during internship at  IsPP BFA INSERM - CNRS UMR 825, under superision of Dr. Leslie Regad & Marine Baillif at Université Paris Cité from Apr to July of 2026. Topic of the internship: Dynamics of HIV-2 protease (PR2)

![License](https://img.shields.io/badge/license-MIT-green) ![Version](https://img.shields.io/badge/version-1.0.0-blue) ![GitHub](https://img.shields.io/badge/GitHub-annhater/M1_STAGE-black?logo=github) ![Build Status](https://img.shields.io/github/actions/workflow/status/annhater/M1_STAGE/ci.yml?branch=main)

## 📋 Table of Contents

- [Features](#features)

## ℹ️ Project Information

- **👤 Author:** annhater
- **📦 Version:** 1.0.0
- **📄 License:** MIT
- **📂 Repository:** [https://github.com/annhater/M1_STAGE](https://github.com/annhater/M1_STAGE)

## Features

The present codes allow calculating the following metrics needed to study PR2 interactions:
  - For active site stability, check: **d<sub>25</sub>** and the codes <code>01_metrics_d25.py</code> or <code>01_metrics_d25.ipynb</code>
  - For flap dynamic, check **d<sub>50</sub>** and the code <code>01_metrics_d25.py</code> or <code>01_metrics_d50.ipynb</code>
  - Other inter-residue distances, such as **d<sub>25-50</sub>** for active site - flap distances can be measured with <code>01_metrics_distance_gen.py</code> 
  - Additionnaly, code to calculate metric of projection **d<sub>50B -> 49A, 50A, 51A</sub>**, to measure the positioning of flap A over B can be found in <code>02_analysis_d50_superimposed.R</code>, written for R
  - Evolution of the conformational dynamic over the course of simulation was measured with **RMSD**. Code: <code>01_metrics_RMSD.py</code>
  - Flexible regions of the structure were defined using **RMSF**. Code: <code>01_metrics_RMSF.py</code> or <code>01_metrics_RMSF.ipynb</code>
  - Change in the phi and psi angles over the course of simulation, calculated with: <code>02_analysis_angles.py</code>
  - To compare similatiory of simulations in different conditions, pairwise RMSD was calculated. Code: <code>02_analysis_pw_interactions.ipynb</code>
  
Additionally, the interactions between the residues of interest were studied:
  - <code>02_analysis_V12_env25_HB.ipynb</code>
  - <code>03_plot_acr_diagram.ipynb</code>
  - <code>02_analysis_env25_heatmap.ipynb</code>
