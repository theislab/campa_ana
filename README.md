# Analysis of 4i data using CAMPA
This repository contains scripts to recreate the models, results, and figures described in 
our preprint 
*"Learning consistent subcellular landmarks to quantify changes in multiplexed protein maps"* 
([Spitzer, Berry et al. (2022]()).
For the main documentation and tutorials, see the official [CAMPA documentation](https://campa.readthedocs.io/).

## Setup
- install and setup [campa](https://github.com/theislab/campa)
- install campa_ana:
  ```
  git clone https://github.com/theislab/campa_ana
  cd campa_ana
  pip install .
  ```
- get data + models
- update `campa.ini` to point to data and experiment folders, 
  and add [NascentRNA](https://github.com/theislab/campa_ana/blob/main/NascentRNA_constants.py) `data_config`.
  The `campa.ini` used to generate this data looked as follows:
  ```
  [DEFAULT]
  data_dir = <path-to-data>
  experiment_dir = <path-to-experiments>

  [data]
  NascentRNA = <path-to-campa_ana>/NascentRNA_constants.py

  [co_occ]
  co_occ_chunk_size = 1e7
  ```

## Data
CSL-derived features from the 184A1 and the HeLa datasets are available [here](https://doi.org/10.6084/m9.figshare.19699651).
Data and pre-trained models are available upon request. 

## Structure of this repository
- `campa_ana`: useful python functions
- `figure_notebooks`: ipython notebooks generating figures
- `figures`: empty folder in which figures will be placed
- `jterator`: param files needed for tissuemaps
- `params`: param files needed for campa
- `R`: useful R functions and R notebooks generating figures
- `scripts`: executable scripts for running campa on HPC
- `workflow`: notebooks explaining how to reproduce datasets, models, and clustering 

## Reproducing results
Notebooks in `workflow` reproduce our results. 
- `01_prepare_data`: download data from TissueMaps and preprocess metadata files
- `02_create_dataset`: create datasets for training and validation of cVAE models
- `03_train`: train cVAE models
- `04_cluster`: cluster latent representation into CSLs and annotate (manually).
- `05_extract_features_all` and `05_extract_features_SBF2`: extract intensity and spatial features from cells using CSLs

For more information on how to use CAMPA and different configuration options for each step, also have a look at the [CAMPA documentation](https://campa.readthedocs.io/).

## Reproducing figures
Notebooks in `figure_notebooks` and `R` reproduce figure panels for the CAMPA paper.

- Figure 1
    - c,d: `R/whole_nucleus_fold_change.R`
    - e,f,g,h: `figure_notebooks/fig1_umap_linear_classifier.ipynb`
- Figure 2
    - a:
    - b,d-i: `R/plot_example_184A1_control_cells.R`
    - c: `R/cluster_loadings_184A1.R`
- Figure 3
    - a,b: `R/fit_mixed_models.R` and `R/make_bubble_plots.R`
    - c-f: `R/plot_example_184A1_meayamycin_cell.R`
    - h,i: `R/plot_co_occurrence_184A1.R`
- Figure 4
    - a-g: `figure_notebooks/fig4_perturbation_comparison.ipynb`
    - h,i: `R/plot_co_occurrence_184A1.R`
    - j: `figure_notebooks/fig4_object_features_CX4561.ipynb`
    - k: `R/plot_example_184A1_cx5461_cell.R`
- Figure 5
    - a,c: `R/plot_example_SBF2_cells.R`
    - b: `R/cluster_loadings_SBF2.R`
    - d: `R/fit_mixed_models_SBF2.R` and `R/make_bubble_plots_SBF2.R`
    - e-h: `R/plot_object_stats_SBF2.R`
- Figure 6
    - b-c: `R/fit_mixed_model_EU_bin.R` and `R/make_bubble_plots_EU_bin.R`
    - d-f: `R/plot_trends_and_examples_EU_heterogeneity.R`

- Supplements:
    - Figure 1: N/A
    - Figure 2 (noise robustness)
        - a: `figure_notebooks/fig1_suppl_noise_robustness.ipynb`
        - b,c: `figure_notebooks/fig1_suppl_cluster_subsampling.ipynb`
    - Figure 3 (perturbation dependent pixel clustering)
        - `figure_notebooks/fig1_umap_linear_classifier.ipynb`
    - Figure 4 (cell-cycle dependent pixel clustering)
        - `figure_notebooks/fig1_suppl_umap_cluster_size_cell_cycle.ipynb`
        - `R/plot_example_184A1_control_cells.R`
    - Figure 5 (cluster loadings)
        - `R/cluster_loadings_184A1.R`
    - Figure 6 (ilastik comparison)
        - `R/compare_CSL_with_ilastik.R`
    - Figure 7 (DMSO bubble-plot)
        - `R/fit_mixed_models_DMSO.R` and `R/make_bubble_plots_DMSO.R`
    - Figure 8 (co-occurrence 184A1)
        - `R/plot_co_occurrence_184A1.R`
    - Figure 9 (co-occurrence 184A1)
        - `R/plot_co_occurrence_184A1.R`
    - Figure 10
    - Figure 11 (cluster loadings SBF2)
        - `R/cluster_loadings_SBF2.R`
    - Figure 12 
        - `R/fit_mixed_models_SBF2.R` and `R/make_bubble_plots_SBF2.R`
        - `R/plot_object_stats_SBF2.R`
    - Figure 13 (object filtering: scatter plots)
        - `figure_notebooks/fig5_object_features_SBF2.ipynb`
    - Figure 14 (object filtering: example cells)
        - `figure_notebooks/fig5_object_features_SBF2.ipynb`
    - Figure 15 (heterogeneity of PML bodies)
        - `figure_notebooks/fig6_pml_umap.ipynb`
