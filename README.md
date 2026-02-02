# Slip-Informed Ground Motion Models in Subduction Zones

![Status](https://img.shields.io/badge/Status-Under%20Review-yellow)
![Language](https://img.shields.io/badge/Language-Python%20%7C%20R-blue)
![License](https://img.shields.io/badge/License-MIT-green)

This repository contains the code, data, and workflows used for the development and analysis of slip-informed Ground Motion Models (GMMs) for subduction zones.

This study proposes new distance metrics ($R_P$, $R_{asp}$) based on Finite Fault Models (FFMs) to better capture seismic source heterogeneity and its impact on ground motion intensity attenuation, comparing them against traditional metrics ($R_{rup}$).

## 📋 Table of Contents
- [Project Description](#project-description)
- [Repository Structure](#repository-structure)
- [Requirements & Installation](#requirements--installation)
- [Workflow & Reproducibility](#workflow--reproducibility)
- [Analyzed Events](#analyzed-events)
- [Data Sources](#data-sources)
- [Citation](#citation)
- [Contact](#contact)

## 📖 Project Description

Traditional Ground Motion Models (GMMs) typically use the closest distance to the rupture ($R_{rup}$) as the primary metric. This work investigates whether incorporating information about the slip distribution within the seismic source improves the prediction of ground motion intensity (PGA, PGV, SA).

Non-linear mixed-effects regressions (`nlmer`) were used to develop both global and regional (Chile, Japan) models. The study analyzes residuals ($dBe$, $dWe$) and compares predictive performance against existing GMMs (e.g., Parker et al., 2020; Montalva et al., 2017).

## 📂 Repository Structure

The repository is organized to ensure reproducibility of the analyses:

```text
├── Data/                     # Databases, slip models, and geographic grids
│   ├── Rp/                   # Pre-calculated Rp and Rp_lock distance values (Excel)
│   ├── Modelos_cosismicos/   # Finite Fault Models (FFM) in .xyz format
│   ├── Modelos_locking/      # Seismic coupling (Locking) grids
│   ├── [Flatfiles CSV]       # Ground motion databases (Drapela et al., PS21)
│   └── [Grid Files]          # Topobathymetry (.nc) and gradients (.int)
│
├── Regressions/              # Regression scripts in R
│   ├── Rrup/                 # Rrup-based models
│   ├── Rp/                   # Rp-based models
│   ├── functions_regression.R # Helper functions for nlmer regression
│   └── Results/              # Resulting coefficients and calculated residuals
│
├── functions/                # Custom Python libraries
│   └── functions_py.py       # Distance calculations, GMM predictions, etc.
│
├── Plots/                    # Python scripts (PyGMT) for figure generation
│   ├── Fig_Residuals_Analysis.py # Statistical analysis of residuals (sigma, tau, phi)
│   ├── Fig_Tohoku_Analysis.py    # Case study analysis: Tohoku 2011
│   ├── Fig_Maule_Analysis.py     # Case study analysis: Maule 2010
│   ├── Fig_Illapel_Analysis.py   # Case study analysis: Illapel 2015
│   └── ...                   # Scripts for other events (Pisagua, Melinka, etc.)
│
└── Figuras_paper/            # Automatic output folder for generated figures (PDF/PNG)
