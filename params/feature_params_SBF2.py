# -- Parameters determining what features are extracted from a given clustering --
feature_params = {
    # -- set up data to extract features from --
    # path to experiment directory relative to EXPERIMENT_DIR
    'experiment_dir': 'VAE_SBF2/CondVAE_siRNA-CC',
    # name of clustering to use
    'cluster_name': None,
    # dir of subsampled clustering to load annotation. Relative to experiment_dir. 
    # Default is taking first of experiment_dir/aggregated/sub-*
    'cluster_dir': "aggregated/sub-0.005_sub-0.33",
    # cluster annotation to use. Defaults to cluster_name
    'cluster_col': None,
    # data dirs to be processed. 
    # Relative to experiment_dir/aggregated/full_data. 
    # If None, all available data_dirs will be processed
    'data_dirs': None,
    #'data_dirs': [  
    #    "HeLa_scrambled/K19", "HeLa_scrambled/L18", "HeLa_scrambled/M19", 
    #    "HeLa_SBF2/K18", "HeLa_SBF2/L19", "HeLa_SBF2/M18"
    #    ],
    # filename to use for saving extracted features.
    'save_name': None,
    # force calculation even when adata exists
    'force': False,
    # -- features to extract --
    # type of features to extract. One or more of intensity, co-occurrence, object-stats
    # Intensity: per-cluster mean and size features. Needs to be calculated first to set up the adata. 
    # Co-occurrence: spatial co-occurrence between pairs of clusters at different distances. 
    # Object stats: number and area of connected components per cluster
    'features': [],
    # parameters for co-occurrence calculation
    'co_occurrence_params': {
        # size of distances interval
        'min': 2.0,
        'max': 679.0,
        'nsteps': 33,
        # use log spacing of co-occurrence intervals
        'logspace': True,
        # number of processes to use to compute co-occurrence scores
        'num_processes': 8
    },
    # parameters for object-stats calculation
    'object_stats_params': {
        # features to extract in mode object-stats
        # possible features: ['area', 'circularity', 'elongation', 'extent']
        'features': ['area', 'circularity', 'elongation', 'extent'],
        # intensity channels to extract mean per cluster from
        'channels': ['20_SP100', '11_PML', '21_COIL', '21_NCL'],
    }
}

# use this list to extract several different features at once.
# final feature params for entry i are obtained by `feature_params.update(variable_feature_params[i])`
variable_feature_params = [
    # intensity features for base clustering seed 1
    {
        'save_name': 'features_seed1.h5ad',
        'cluster_name': 'clustering_res0.9_sub-0.33_seed1',
        'features': ['intensity'],
    },
    # intensity features for based clustering seed 2
    {
        'save_name': 'features_seed2.h5ad',
        'cluster_name': 'clustering_res0.9_sub-0.33_seed2',
        'features': ['intensity'],
    },
    # intensity features for based clustering seed 3
    {
        'save_name': 'features_seed3.h5ad',
        'cluster_name': 'clustering_res0.9_sub-0.33_seed3',
        'features': ['intensity'],
    },
    # intensity + co-occurrence + object stats for annotated clustering based on seed 3
    {
        'save_name': 'features_seed3_annotation.h5ad',
        'cluster_name': 'clustering_res0.9_sub-0.33_seed3',
        'cluster_col': 'annotation',
        'features': ['intensity', 'co-occurrence', 'object-stats'],
    },
    # intensity + co-occurrence for nucleus/cytoplasm grouping for seed 3
    {
        'save_name': 'features_seed3_annotation_cytoplasm.h5ad',
        'cluster_name': 'clustering_res0.9_sub-0.33_seed3',
        'cluster_col': 'annotation_cytoplasm',
        'features': ['intensity', 'co-occurrence'],
    },
    # -- old feature params for reproducibility below --
    # NOTE: initially co-occ scores were calculated for distances 2-320, and 320-679 separately.
    # when re-running, just execute the above calls to directly calculate all distances.
    # below are the calls for calculating distances individually. They were merged using code in workflow/05_extract_features.ipynb. 

    # end-distance chosen such that spacing between samples in log-space is similar to spacing before. 679 ~= 2**(np.log2(320)+ 0.27118252*4)
    # spacing before: np.log2(np.logspace(np.log2(2),np.log2(320),28, base=2))[1:] - np.log2(np.logspace(np.log2(2),np.log2(320),28, base=2))[:-1]
    # spacing now: samples = np.concatenate([np.logspace(np.log2(2),np.log2(320),28, base=2).astype(np.float32), np.logspace(np.log2(320),np.log2(679),5, base=2).astype(np.float32)])
    # np.log2(samples)[1:] - np.log2(samples)[:-1]
    # co-occurrence for annotation for seed 3 (small distances)
    #{
    #    'save_name': 'features_seed3_annotation.h5ad',
    #    'cluster_name': 'clustering_res0.9_sub-0.33_seed3',
    #    'cluster_col': 'annotation',
    #    'features': ['intensity', 'co-occurrence'],
    #    'co_occurrence_params': {
    #        'max': 320,
    #        'nsteps': 28,
    #    },
    #},
    # co-occurrence for annotation for seed 3 (large distances)
    #{
    #    'save_name': 'features_seed3_annotation2.h5ad',
    #    'cluster_name': 'clustering_res0.9_sub-0.33_seed3',
    #    'cluster_col': 'annotation',
    #    'features': ['intensity', 'co-occurrence'],
    #    'co_occurrence_params': {
    #        'min': 320,
    #        'max': 679,
    #        'nsteps': 5
    #    },
    #},
    # co-occurrence for nucleus/cytoplasm grouping for seed 3 (small distances)
    #{
    #    'save_name': 'features_seed3_annotation_cytoplasm.h5ad',
    #    'cluster_name': 'clustering_res0.9_sub-0.33_seed3',
    #    'cluster_col': 'annotation_cytoplasm',
    #    'features': ['intensity', 'co-occurrence'],
    #    'co_occurrence_params': {
    #        'max': 320,
    #        'nsteps': 28,
    #    },
    #},
    # co-occurrence for nucleus/cytoplasm grouping for seed 3 (large distances)
    #{
    #    'save_name': 'features_seed3_annotation_cytoplasm2.h5ad',
    #    'cluster_name': 'clustering_res0.9_sub-0.33_seed3',
    #    'cluster_col': 'annotation_cytoplasm',
    #    'features': ['intensity', 'co-occurrence'],
    #    'co_occurrence_params': {
    #        'min': 320,
    #        'max': 679,
    #        'nsteps': 5
    #    },
    #},
]