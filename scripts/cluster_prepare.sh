#! /bin/bash

#SBATCH -o cluster_prepare_%j.out
#SBATCH -e cluster_prepare_%j.out
#SBATCH -J cluster_prepare
#SBATCH -p cpu_p
#SBATCH -c 1
#SBATCH --mem 255G
##SBATCH --exclude=ibis-ceph-0[02-05,08-19],ibis216-010-0[20-37,64],icb-rsrv[05-06,08]
#SBATCH -t 20:00:00
#SBATCH --nice=10000

source ~/.bashrc
conda activate campa

# VAE_all experiment
# create subsampled mpp_cluster data
#for exp in VAE CondVAE_pert-CC MPPleiden; do
#    campa cluster VAE_all/$exp create --subsample --frac 0.001 --save-dir aggregated/sub-0.001
#done

# prepare full dataset for projecting clustering to
#for exp in VAE CondVAE_pert-CC MPPleiden; do
#    campa cluster VAE_all/$exp prepare-full --save-dir aggregated/full_data
#done

# VAE_SBF2 experiment
# create subsampled mpp_cluster data
#for exp in VAE CondVAE_siRNA-CC MPPleiden; do
#    campa cluster VAE_SBF2/$exp create --subsample --frac 0.005 --save-dir aggregated/sub-0.005
#done

# prepare full dataset for projecting clustering to
#for exp in VAE CondVAE_siRNA-CC MPPleiden; do
#    campa cluster VAE_SBF2/$exp prepare-full --save-dir aggregated/full_data
#done

# ablation studies
for exp in CondVAE_pert-CC_noneigh CondVAE_pert-CC_neigh5 CondVAE_pert-CC_neigh7; do
    campa cluster VAE_all/$exp create --subsample --frac 0.001 --save-dir aggregated/sub-0.001
    # TODO maybe not needed, might be possible to manually predict this for the few example cells needed
    #campa cluster VAE_all/$exp prepare-full --save-dir aggregated/full_data
done