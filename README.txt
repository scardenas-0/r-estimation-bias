R Estimation Bias Figures
=========================

This repository contains the MATLAB code used to generate the main-text figures and the Figure S.2 sample-size validation analysis.

Main-text and supplemental figures:
1. Open MATLAB and set this folder as the current folder.
2. Run:

   generate_data
   make_figures

`generate_data.m` creates the MATLAB data files in `data/`, and `make_figures.m` reads those files and exports the figures to `figures/`.

Figure S.2 sample-size validation:
1. Open MATLAB and set this folder as the current folder.
2. Run:

   generate_figure_S2_sample_size_validation

Required MATLAB toolbox:
Statistics and Machine Learning Toolbox

The `figures` folder contains the exported figures.

The `data` folder contains CSV summaries and MATLAB result files for the Figure S.2 analysis.
