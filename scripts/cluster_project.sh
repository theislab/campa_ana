#! /bin/bash

#SBATCH -o cluster_project_%j.out
#SBATCH -e cluster_project_%j.out
#SBATCH -J cluster_project
#SBATCH -p cpu_p
#SBATCH -c 1
#SBATCH --mem 255G
##SBATCH --exclude=ibis-ceph-0[02-05,08-19],ibis216-010-0[20-37,64],icb-rsrv[05-06,08]
#SBATCH -t 20:00:00
#SBATCH --nice=10000

source ~/.bashrc
conda activate pelkmans-3.9

BASE=/home/icb/hannah.spitzer/projects/pelkmans/software_new
SCRIPT=$BASE/campa/cli/cluster.py

# VAE_all experiment
# project clustering
python $SCRIPT VAE_all/CondVAE_pert-CC project aggregated/sub-0.001 --save-dir aggregated/full_data --cluster-name "clustering_res0.5"


# VAE_SBF2 experiment
# project clustering
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC project aggregated/sub-0.005_sub-0.33 --save-dir aggregated/full_data --cluster-name "clustering_res0.9_sub-0.33_seed1"
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC project aggregated/sub-0.005_sub-0.33 --save-dir aggregated/full_data --cluster-name "clustering_res0.9_sub-0.33_seed2"
#python $SCRIPT VAE_SBF2/CondVAE_siRNA-CC project aggregated/sub-0.005_sub-0.33 --save-dir aggregated/full_data --cluster-name "clustering_res0.9_sub-0.33_seed3"
