#! /bin/bash

#SBATCH -o train_%j.out
#SBATCH -e train_%j.out
#SBATCH -J train
#SBATCH -p gpu_p
#SBATCH --gres=gpu:1
#SBATCH -c 1
#SBATCH --mem 64G
#SBATCH --qos=gpu
#SBATCH -t 12:00:00
#SBATCH --nice=10000
#SBATCH --exclude=icb-gpusrv0[1-2]

# set up
source ~/.bashrc
conda activate campa

# run script
BASE=/home/icb/hannah.spitzer/projects/pelkmans/software_new

# train, evaluate, and compare all experiments
campa train trainval --config $BASE/campa_ana/params/experiment_params_all.py
#campa train all --config $BASE/campa_ana/params/experiment_params_SBF2.py

# code for evaluation + comparison only:
#campa train evaluate --experiment-dir VAE_all
campa train compare --experiment-dir VAE_all # --exp-name MPPleiden VAE CondVAE_pert-CC 
#campa train evaluate --experiment-dir VAE_SBF2
#campa train compare --experiment-dir VAE_SBF2

