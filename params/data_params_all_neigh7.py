# --- Parameters for creating dataset for all perturbations ---
# conditions: perturbation and cell cycle
import os

data_params = {
    # name of the resulting dataset that is defined by these params (relative to DATA_DIR/datasets)
    'dataset_name': '184A1_all_frac0005_neigh7_cond_pert-CC',
    # name of data config (registered in config.ini)
    'data_config': "NascentRNA",
    # where to read data from (relative to DATA_DIR defined in data config)
    'data_dirs': 
            [os.path.join('184A1_unperturbed', well) for well in ['I09', 'I11', 'J10', 'J12']] + \
            [os.path.join('184A1_DMSO', well) for well in ['I14', 'J16']] + \
            [os.path.join('184A1_AZD4573', well) for well in ['I13', 'I17', 'J14', 'J18', 'J21']] + \
            [os.path.join('184A1_CX5461', well) for well in ['I18', 'J09', 'J22']] + \
            [os.path.join('184A1_TSA', well) for well in ['I16', 'J13', 'J20']] + \
            [os.path.join('184A1_triptolide', well) for well in ['I10', 'J15']] + \
            [os.path.join('184A1_meayamycin', well) for well in ['I12', 'I20']],
    'channels': [
        '01_CDK9_pT186', '01_PABPC1', '02_CDK7', '03_CDK9', '03_RPS6', '05_GTF2B', '05_Sm', '07_POLR2A', '07_SETD1A', 
        '08_H3K4me3', '09_CCNT1', '09_SRRM2', '10_H3K27ac', '10_POL2RA_pS2', '11_KPNA2_MAX', '11_PML', '12_RB1_pS807_S811', 
        '12_YAP1', '13_PABPN1', '13_POL2RA_pS5', '14_PCNA', '15_SON', '15_U2SNRNPB', '16_H3', '17_HDAC3', '17_SRSF2', 
        '18_NONO', '19_KPNA1_MAX', '20_ALYREF', '20_SP100', '21_COIL', '21_NCL', '00_DAPI', '07_H2B'
    ], 
    # list of conditions. Should be defined in data config. 
    # The suffix '_one_hot' will convert the condition in a one-hot encoded vector.
    # Conditions are concatenated, except when they are defined as a list of lists. 
    # In this case the condition is defined as a pairwise combination of the conditions.
    'condition': ['perturbation_duration_one_hot', 'cell_cycle_one_hot'],
    'condition_kwargs': {
        'cond_params': {}
    },
    # train/val/test split
    'split_kwargs': {
        'train_frac': 0.8,
        'val_frac': 0.1,
    },
    'test_img_size': 225,
    # subset to objects with certain metadata.
    'subset': True,
    # kwargs to MPPData.subset() defining which object to subset to
    'subset_kwargs': {
        'frac': None, # special kwarg for random subsetting of objects
        'nona_condition': True,  # special kwarg for removing all objects with NAN condition
        'cell_cycle': 'NO_NAN'
    },
    # subsampling of pixels (only for train/val)
    'subsample': True,
    # kwargs for MPPData.subsample() defining the fraction of pixels to be sampled
    'subsample_kwargs': {
        'frac': None,
        'frac_per_obj': 0.005,
        'num': None,
        'num_per_obj': None,
    },
    # neighborhood information
    'neighborhood': True,
    'neighborhood_size': 7,
    # normalisation
    'normalise': True,
    'normalise_kwargs': {
        # background_value is column name in CHANNELS_METADATA, or list of floats per channel
        'background_value': 'mean_background',
        'percentile': 98.0,
        'rescale_values': [],
    },
    # make results reproducible
    'seed': 42,
}
