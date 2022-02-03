#! /bin/bash

#SBATCH -o features_%j.out
#SBATCH -e features_%j.out
#SBATCH -J features
#SBATCH -p cpu_p
#SBATCH -c 8
#SBATCH --mem 127G
#SBATCH -t 40:00:00
#SBATCH --nice=10000

source ~/.bashrc
conda activate pelkmans-3.9

BASE=/home/icb/hannah.spitzer/projects/pelkmans/software_new
SCRIPT=$BASE/campa/cli/extract_features.py


# -- VAE_SBF2 experiment --
# extract intensity features for each clustering
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed1" --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed1.h5ad --force intensity
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed1" --cluster-dir aggregated/sub-0.005_sub-0.33 --cluster-col annotation --save-name features_seed1_annotation.h5ad --force intensity
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed2" --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed2.h5ad --force intensity
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3.h5ad --force intensity
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" --cluster-dir aggregated/sub-0.005_sub-0.33 --cluster-col annotation --save-name features_seed3_annotation.h5ad --force intensity

# object stats (for annotated seed3)
python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" --cluster-col annotation --save-name features_seed3_annotation.h5ad object-stats --area-threshold 10

# NOTE: co-occurrence calculation for this experiment takes about 40 hours per well
# co-occurrence (for annotated seed3)
#python $SCRIPT \
#    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" co-occurrence \
#    --cluster-col annotation --co-logspace --co-nsteps 33 \
#    --co-minval 2 --co-maxval 679  --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation.h5ad \
#    --num-processes 8 \
#    --data-dir "HeLa_scrambled/K19" "HeLa_scrambled/L18" "HeLa_scrambled/M19" "HeLa_SBF2/K18" "HeLa_SBF2/L19" "HeLa_SBF2/M18"
    
# co-occurrence (for nucleus/cytoplasm grouping of seed3)
#python $SCRIPT \
#    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" co-occurrence \
#    --cluster-col "annotation_cytoplasm" --co-logspace --co-nsteps 33 \
#    --co-minval 2 --co-maxval 679  --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation_cytoplasm.h5ad \
#    --num-processes 8 \
#    --data-dir "HeLa_scrambled/K19" "HeLa_scrambled/L18" "HeLa_scrambled/M19" "HeLa_SBF2/K18" "HeLa_SBF2/L19" "HeLa_SBF2/M18"

# NOTE: initially co-occ scores were calculated for distances 2-320, and 320-679 separately.
# when re-running, just execute the above calls to directly calculate all distances.
# below are the calls for calculating distances individually. They were merged using code in workflow/05_extract_features.ipynb. 

# end-distance chosen such that spacing between samples in log-space is similar to spacing before. 679 ~= 2**(np.log2(320)+ 0.27118252*4)
# spacing before: np.log2(np.logspace(np.log2(2),np.log2(320),28, base=2))[1:] - np.log2(np.logspace(np.log2(2),np.log2(320),28, base=2))[:-1]
# spacing now: samples = np.concatenate([np.logspace(np.log2(2),np.log2(320),28, base=2).astype(np.float32), np.logspace(np.log2(320),np.log2(679),5, base=2).astype(np.float32)])
# np.log2(samples)[1:] - np.log2(samples)[:-1]

# small distances for annotated seed3
#python $SCRIPT \
#    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" co-occurrence \
#    --cluster-col annotation --co-logspace --co-nsteps 28 \
#    --co-maxval 320  --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation.h5ad \
#    --num-processes 8 \
#    --data-dir "HeLa_scrambled/K19" "HeLa_scrambled/L18" "HeLa_scrambled/M19" "HeLa_SBF2/K18" "HeLa_SBF2/L19" "HeLa_SBF2/M18"

# large distances for annotated seed3
#python $SCRIPT \
#    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" co-occurrence \
#    --cluster-col annotation --co-logspace --co-nsteps 5 \
#    --co-minval 320 --co-maxval 679 --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation2.h5ad \
#    --num-processes 8 \
#    --data-dir "HeLa_scrambled/K19" "HeLa_scrambled/L18" "HeLa_scrambled/M19" "HeLa_SBF2/K18" "HeLa_SBF2/L19" "HeLa_SBF2/M18"

# small distances for for nucleus/cytoplasm grouping of seed3
#python $SCRIPT \
#    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" co-occurrence \
#    --cluster-col "annotation_cytoplasm" --co-logspace --co-nsteps 28 \
#    --co-maxval 320  --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation_cytoplasm.h5ad \
#    --num-processes 8 \
#    --data-dir "HeLa_scrambled/K19" "HeLa_scrambled/L18" "HeLa_scrambled/M19" "HeLa_SBF2/K18" "HeLa_SBF2/L19" "HeLa_SBF2/M18"

# large distances for for nucleus/cytoplasm grouping of seed3
#    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" co-occurrence \
#    --cluster-col "annotation_cytoplasm" --co-logspace --co-nsteps 5 \
#    --co-minval 320 --co-maxval 679 --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation_cytoplasm2.h5ad \
#    --num-processes 8 \
#    --data-dir "HeLa_scrambled/K19" "HeLa_scrambled/L18" "HeLa_scrambled/M19" "HeLa_SBF2/K18" "HeLa_SBF2/L19" "HeLa_SBF2/M18"