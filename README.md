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

- **theta2**: creates the theta2-plot, e.g. when changing thresholds of `gamma_prediction` or the theta2-cut
- **skymap**: creates a 2-dimensional skymap of the reconstructed source position via `skymaps.py`
- **protonskymap**: does the above for the proton MC simulations (models have to be applied to the simulations first)
- **gammaskymap**: does the above for the gamma MC simulations (models have to be applied to the simulations first)
- **theta**: does both of the above
- **mcdata\_p**: uses `plot_data_mc_comparison.py` to generate comparison plots of the features of observed data and the proton MC simulations.
- **mcdata**: uses `ks_test.py` to generate data MC comparison plots for observed data, gamma simulations and proton simulations.
- **plot\_performance**: uses the aict-tools command line functions to generate performance plots for the machine learning models (gamma hadron separation and energy reconstruction) generated on the simulations.
- **meta**: adds the meta data of the open data sample observations to the runs group of the hdf5 files via `add_meta_data.py`. Needs to be applied before starting the analysis.
- **rate**: calls `plot_ratescan.py` generate a plot of the number of observation events as a function of a cut on the parameter 'size'
- **delta**: calls `delta_delta/plot_delta_delta.py` to generate a plot showing the difference of the reconstructed `delta` and the expected `delta` from the known source position to investigate origin reconstruction.
- **sdelta**: same as the above, but calling `delta_delta/plot_specific_delta_delta.py` to compare the FACT-Tools `delta` with the exact same events of the PhotonStream.
- **cuts**: uses the [performance python analysis'](https://github.com/fact-project/performance-paper-analysis) `detection_gridsearch.py` to calculate the optimal `gamma\_prediction` and theta2-cut threshold for the maximum detection significance

## Additional Investigations

There are a few additional investigations performed to understand the PhotonStream data better. They
are basically represented by the different directories in this repository and will be explained in
the following. The investigations generally consist of a dedicated Makefile, which produces a number
of plots within their respective directory.

### camera\_images\_cleaning

The DBSCAN cleaning within the PhotonStream is a new way to clean IACT images within the
three-dimensional point cloud. Therefore, it is interesting to investigate the differences.
`plot_cleaned_images.py` plots the camera images and, if desired, the cleaned pixels within the
two-dimensional camera plane for the classical FACT-Tools cleaning on both data representations
plus the DBSCAN cleaning on PhotonStream data. Further tunings and information can be found in
the python script.

### delta\_delta

To further understand the poor out\-of\-the\-box performance of the origin reconstruction on
PhotonStream data, the feature `delta` can be investigated. It represents the orientation of the
air-shower within the camera with respect to the x-axis. Since the source position of the Crab
Nebula is well known and with the help of the telescope's pointing position can be reconstructed
within the camera, the difference between the calculated `delta` and this expected `delta` can be
calculated. This can be done with `plot_delta_delta.py`. Furthermore, the difference between the
FACT-Tools `delta` and the PhotonStream `delta` can be calculated by `plot_specific_delta_delta.py`.
This way, the difference of this new data representation to the classical Largest Pulse data can be
examined.

### eps\_comparison\_plots

To investigate the effects of the new cleaning and whether it can be optimized, two parameters can
be varied. The minimum number of photons within a dense cluster (m) and the minimum distance of two
photons to be considered dense (epsilon). In this work, the parameter epsilon can be investigated.
Especially with respect to its impact on mismatches between data and simulations. To do so the
mismatches can be investigated directly for samples generated with different epsilon. Furthermore,
via `plot_eps_scores.py` the high level performance scores of the models and the analysis results
can b einvestigated as a function of the epsilon. Unfortunately, there is no complete automation of
this, since the performance results are only printed to the terminal output and therefore are
hard-coded here.

### pe\_difference

By using the new single photon extraction the resulting images of the PhotonStream differ from the
classical representation. To investigate the differences in the number of reconstructed
photon-equivalents, this Makefile can be used. `plot_cleaned_images.py` can be used to plot the
two-dimensional camera images for the two data representations available in the open data sample.
Additionally to those two images, the difference in PE is plotted into the camera image in the middle.
Furthermore, `plot_histograms.py` creates histograms for the PE differences per pixel and the
differences of the mean PE per event.

### time\_slices

As a first investigation of the new time features, the PhotonStream is offering, the distributions of
the arrival times of single photons can be displayed via `slices.py`. Furthermore, since only bright
pixels are interesting as signal, a cut on the number of minimum PE can be set and only such pixel
arrival times put in the histograms for observations as well as simulations.
