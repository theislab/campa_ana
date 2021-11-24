#! /bin/bash

#SBATCH -o features_%j.out
#SBATCH -e features_%j.out
#SBATCH -J features
#SBATCH -p cpu_p
#SBATCH -c 1
#SBATCH --mem 255G
##SBATCH --exclude=ibis-ceph-0[02-05,08-19],ibis216-010-0[20-37,64],icb-rsrv[05-06,08]
#SBATCH -t 10:00:00
#SBATCH --nice=10000

source ~/.bashrc
conda activate pelkmans-3.9

BASE=/home/icb/hannah.spitzer/projects/pelkmans/software_new
SCRIPT=$BASE/miann/cli/extract_features.py

# VAE_all experiment
#python $SCRIPT 

# VAE_SBF2 experiment
# extract intensity features for each clustering
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed1" --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed1.h5ad --force intensity
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed1" --cluster-dir aggregated/sub-0.005_sub-0.33 --cluster-col annotation --save-name features_seed1_annotation.h5ad --force intensity
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed2" --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed2.h5ad --force intensity
python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC "clustering_res0.9_sub-0.33_seed3" --cluster-dir aggregated/sub-0.005_sub-0.33 --save-name features_seed3.h5ad --force intensity
