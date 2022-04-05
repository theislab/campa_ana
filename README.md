# Analysis of 4i data using CAMPA
This repository contains scripts to recreate the models, results, and figures described in the CAMPA paper.

## Setup
- install and setup [campa](https://github.com/theislab/campa)
- install campa_ana:
```
git clone https://github.com/theislab/campa_ana
cd campa_ana
pip install .
```
- download data + models (TODO)
- update `camap.ini` to point to data and experiment folders

## Data
TODO - how should data be provided and accessed. Which data to provide?

## Structure of this repository
- `campa_ana`: useful python functions
- `figure_notebooks`: ipython notebooks generating figures
- `figure_notebooks_R`: R notebooks generating figures
- `figures`: empty folder in which figures will be placed
- `jterator`: param files needed for tissuemaps
- `params`: param files needed for campa
- `R`: useful R functions
- `scripts`: executable scripts for running campa on HPC
- `workflow`: notebooks explaining how to reproduce datasets, models, and clustering 

## Reproducing results
Notebooks in `workflow` reproduce our results. 
- `01_prepare_data`: download data from TissueMaps and preprocess metadata files
- `02_create_dataset`: create datasets for training and validation of cVAE models
- `03_train`: train cVAE models
- `04_cluster`: cluster latent representation into CSLs and annotate (manually).
- `05_extract_features_all` and `05_extract_features_SBF2`: extract intensity and spatial features from cells using CSLs

For more information on how to use CAMPA and different configuration options for each step, also have a look at the [CAMPA documentation](TODO).

## Reproducing figures
Notebooks in `figure_notebooks` reproduce figure panels for the CAMPA paper.

TODO: add info on which notebooks output which panels

- Figure 1
    - C,D: Scott
    - E,F,G,H: `fig1_umap_linear_classifier.ipynb`

- Supplements:
    - Figure 2 (noise robustness)
        - A: `fig1_suppl_noise_robustness.ipynb`
        - B,C: `fig1_suppl_cluster_subsampling.ipynb`
    - Figure 3 (perturbation dependent pixel clustering)
        - `fig1_umap_linear_classifier.ipynb`
        
## Figure files
- figures can be found [in this nextcloud folder](https://hmgubox2.helmholtz-muenchen.de/index.php/s/36PrZwt4cMniLfW)
- sync figures: `rsync -rva hpc-submit01.scidom.de:/home/icb/hannah.spitzer/projects/pelkmans/software_new/campa_ana/figures/ ~/projects/pelkmans/software_new/campa_ana/figures --exclude=".*"`
