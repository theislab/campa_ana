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
SCRIPT=$BASE/miann/cli/extract_features.py


# -- VAE_SBF2 experiment --
# extract intensity features for each clustering
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed1" --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed1.h5ad --force intensity
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed1" --cluster-dir aggregated/sub-0.005_sub-0.33 --cluster-col annotation --save-name features_seed1_annotation.h5ad --force intensity
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed2" --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed2.h5ad --force intensity
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3.h5ad --force intensity
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" --cluster-dir aggregated/sub-0.005_sub-0.33 --cluster-col annotation --save-name features_seed3_annotation.h5ad --force intensity

# co-occurrence (for annotated seed3)
#python $SCRIPT \
#    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" co-occurrence \
#    --cluster-col annotation --co-logspace --co-nsteps 28 \
#    --co-maxval 320  --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation.h5ad \
#    --num-processes 8 \
#    --data-dir "HeLa_scrambled/K19"  
    
#python $SCRIPT \
#    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" intensity co-occurrence \
#    --cluster-col "annotation_cell" --co-logspace --co-nsteps 28 \
#    --co-maxval 320  --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation_cell.h5ad \
#    --num-processes 8 \
#    --data-dir "HeLa_SBF2/K18" "HeLa_SBF2/L19" "HeLa_SBF2/M18"
    
python $SCRIPT \
    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" intensity co-occurrence \
    --cluster-col "annotation_cytoplasm" --co-logspace --co-nsteps 28 \
    --co-maxval 320  --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation_cytoplasm.h5ad \
    --num-processes 8 \
    --data-dir "HeLa_scrambled/K19" "HeLa_scrambled/L18" "HeLa_scrambled/M19"
    
    #"HeLa_scrambled/K19" "HeLa_scrambled/L18" "HeLa_scrambled/M19"

    #"HeLa_SBF2/K18" "HeLa_SBF2/L19" "HeLa_SBF2/M18"
    
# K18 running - 4811021 - restarted - 4824859
# L19 running - 4813010 - restarted - 4824868
# M18 running - 4818623 - DONE
# K19 running - 4818649 - restarted - 4830351
# L18 running - 4818854 - restarted - 4830308
# M19 running - 4818855 - restarted - 4824887
# was looking at 05_extract features to judge if co-occ is finished

#    "HeLa_scrambled/K19" "HeLa_scrambled/L18" "HeLa_scrambled/M19"

    #"HeLa_SBF2/K18" "HeLa_SBF2/L19" "HeLa_SBF2/M18"
