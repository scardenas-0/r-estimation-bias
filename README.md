# Bias in Estimating Subcritical Reproduction Numbers Under Imperfect Observation and Overlapping Transmission Chains: Theoretical Framework and Application to Mpox

Seth Blumberg, Santiago D. Cardenas, Andrew W. Liu, Taye Samuel Faniran, Juliet R. C. Pulliam, James O. Lloyd-Smith

## Repository information

This repository accompanies the article "Bias in Estimating Subcritical Reproduction Numbers Under Imperfect Observation and Overlapping Transmission Chains: Theoretical Framework and Application to Mpox" (Blumberg et al. 2026)

This repository contains the MATLAB code used to generate the main-text figures and the Figure S.2 sample-size validation analysis.

## Main-text and supplemental figures

To regenerate the main-text and supplement and figures, run the scripts in this order:

```matlab
generate_data
make_figures
generate_figure_S2_sample_size_validation
```
`generate_data.m` creates the MATLAB data files in `data/`, and `make_figures.m` reads those files and exports the figures to `figures/`. `generate_figure_S2_sample_size_validation` creates Figure S.2.

## Requirement

MATLAB Statistics and Machine Learning Toolbox.

## Outputs

The `figures` folder contains the exported figures.

The `data` folder contains CSV summaries and MATLAB result files for the Figure S.2 analysis.
