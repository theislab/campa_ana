# Analysis of 4i data using CAMPA
This repository contains scripts to recreate the models, results, and figures described in the 4i paper.

We provide both data and pre-trained models described in the paper. 

## Setup
- install campa + campa_ana
- download data + models
- change `config.ini` to point to data

## Data
TODO - how should data be provided and accessed. Which data to provide?

## Workflow
Notebooks in `workflow` reproduce our results. 
- `01_prepare_data`: download data from TissueMaps and preprocess metadata files
- `02_create_dataset`: create datasets for training and validation of cVAE models
- `03_train`: train cVAE models
- `04_cluster`: cluster latent representation into CSLs and annotate (manually)
- `05_extract_features_all` and `05_extract_features_SBF2`: extract intensity and spatial features from cells using CSLs

For more information on how to use CAMPA and different configuration options for each step, also have a look at the tutorials provided together with CAMPA (LINK)