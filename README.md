# Open Data PhotonStream Analysis
Analysis of the FACT Open Data Sample's PhotonStream Data.

This repository contains all scripts used for my
[Master thesis](https://github.com/KevSed/Master-Thesis) on the new PhotonStream data representation.
The analysis is based on the classical FACT analysis and adapted to the existing
[FACT-Tools](https://github.com/fact-project/fact-tools) analysis of the FACT open
[data](https://fact-project.org/data/) crab sample, which can be found
[here](https://github.com/fact-project/open_crab_sample_analysis/).

The PhotonStream data representation is described in detail
[here](https://github.com/fact-project/photon_stream). It contains the photon extraction, calibration
and cleaning of the data. The generation of the features for the analysis is implemented in the
[FeatureStream](https://github.com/KevSed/FeatureStream), alongside with additional cleanings and
new time features.

## Main Analysis Tasks
All the main analysis tasks are performed in the right order by the `Makefile` in the root directory
of this repository. It will start by generating features from the data files available on the FACT
open data website. The previous download of these files is therefore necessary. The python script
`gendata.py` will start creating the features using multiprocessing. It is currently set to use 48
cores. It furthermore contains the possibility to change this number and the clustering parameter
of the cleaning via command line arguments.

The application of quality cuts as defined in `configs/quality_cuts.yaml`, the splitting of the data
samples into training and test samples, the training of the random forest models on the simulations,
the application of these models to the data (and if desired the simulations), and the generation of
the resulting theta2-plot are performed by the great artificial intelligence tools for Cherenkov
telescopes, [aict-tools](https://github.com/fact-project/aict-tools).

For further investigations of the PhotonStream, a couple of additional plots can be generated. The
following scripts can be called via the Makefile.

- theta2: creates the theta2-plot, e.g. when changing thresholds of `gamma_prediction` or the theta2-cut
- skymap: creates a 2-dimensional skymap of the reconstructed source position via `skymaps.py`
- protonskymap: does the above for the proton MC simulations (models have to be applied to the simulations first)
- gammaskymap: does the above for the gamma MC simulations (models have to be applied to the simulations first)
- theta: does both of the above
- mcdata\_p: uses `plot_data_mc_comparison.py` to generate comparison plots of the features of observed data and the proton MC simulations.
- mcdata: uses `ks_test.py` to generate data MC comparison plots for observed data, gamma simulations and proton simulations.
- plot\_performance: uses the aict-tools command line functions to generate performance plots for the machine learning models (gamma hadron separation and energy reconstruction) generated on the simulations.
- meta: adds the meta data of the open data sample observations to the runs group of the hdf5 files via `add_meta_data.py`. Needs to be applied before starting the analysis.
- rate: calls `plot_ratescan.py` generate a plot of the number of observation events as a function of a cut on the parameter 'size'
- delta: calls `delta_delta/plot_delta_delta.py` to generate a plot showing the difference of the reconstructed `delta` and the expected `delta` from the known source position to investigate origin reconstruction.
- sdelta: same as the above, but calling `delta_delta/plot_specific_delta_delta.py` to compare the FACT-Tools `delta` with the exact same events of the PhotonStream.
- cuts: uses the [performance python analysis'](https://github.com/fact-project/performance-paper-analysis) `detection_gridsearch.py`
