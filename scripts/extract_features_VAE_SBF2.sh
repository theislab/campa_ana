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

# NOTE: co-occurrence calculation for this experiment takes about 40 hours per well
# co-occurrence (for annotated seed3)
#python $SCRIPT \
#    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" co-occurrence \
#    --cluster-col annotation --co-logspace --co-nsteps 28 \
#    --co-maxval 320  --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation.h5ad \
#    --num-processes 8 \
#    --data-dir "HeLa_scrambled/K19" "HeLa_scrambled/L18" "HeLa_scrambled/M19" "HeLa_SBF2/K18" "HeLa_SBF2/L19" "HeLa_SBF2/M18"
    
# co-occurrence (for cell cluster)
# NOTE: extracting cell only results in co-occ of 1 everywhere (kinda makes sense, need another cluster to compare to)
#python $SCRIPT \
#    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" co-occurrence \
#    --cluster-col "annotation_cell" --co-logspace --co-nsteps 28 \
#    --co-maxval 320  --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation_cell.h5ad \
#    --num-processes 8 \
#    --data-dir "HeLa_scrambled/K19" "HeLa_scrambled/L18" "HeLa_scrambled/M19" "HeLa_SBF2/K18" "HeLa_SBF2/L19" "HeLa_SBF2/M18"
    
# co-occurrence (for nucleus/cytoplasm grouping of seed3)
python $SCRIPT \
    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" co-occurrence \
    --cluster-col "annotation_cytoplasm" --co-logspace --co-nsteps 28 \
    --co-maxval 320  --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation_cytoplasm.h5ad \
    --num-processes 8 \
    --data-dir "HeLa_SBF2/M18"

# NOTE calculate co-occ scores for larger distances. save as separate file, need to merge files later. 
# end-distance chosen such that spacing between sammples in log-space is similar to spacing before. 679 ~= 2**(np.log2(320)+ 0.27118252*4)
# spacing before: np.log2(np.logspace(np.log2(2),np.log2(320),28, base=2))[1:] - np.log2(np.logspace(np.log2(2),np.log2(320),28, base=2))[:-1]
# spacing now: samples = np.concatenate([np.logspace(np.log2(2),np.log2(320),28, base=2).astype(np.float32), np.logspace(np.log2(320),np.log2(679),5, base=2).astype(np.float32)])
# np.log2(samples)[1:] - np.log2(samples)[:-1]

#python $SCRIPT \
#    VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" intensity co-occurrence \
#    --cluster-col "annotation_cytoplasm" --co-logspace --co-nsteps 5 \
#    --co-minval 320 --co-maxval 679  --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3_annotation_cytoplasm2.h5ad \
#    --num-processes 8 \
#    --data-dir "HeLa_SBF2/L19"

    # running: M18, L19


