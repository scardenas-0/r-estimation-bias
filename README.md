# R Estimation Bias Figures

This repository contains the MATLAB code used to generate the main-text figures and the Figure S.2 sample-size validation analysis.

## Main-text and supplemental figures

To regenerate the main-text and supplement (except Figure S.2) results and figures, run the scripts in this order:

```matlab
generate_data
make_figures
```

`generate_data.m` creates the MATLAB data files in `data/`, and `make_figures.m` reads those files and exports the figures to `figures/`.

## Figure S.2 sample-size validation

To regenerate Figure S.2, open this folder in MATLAB and run:

```matlab
generate_figure_S2_sample_size_validation
```

## Requirement

MATLAB Statistics and Machine Learning Toolbox.

## Outputs

The `figures` folder contains the exported figures.

The `data` folder contains CSV summaries and MATLAB result files for the Figure S.2 analysis.
