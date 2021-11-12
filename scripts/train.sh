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
conda activate pelkmans-3.9

# run script
BASE=/home/icb/hannah.spitzer/projects/pelkmans/software_new
SCRIPT=$BASE/miann/cli/train.py

#python $SCRIPT all --config $BASE/miann_ana/params/experiment_params_all.py
python $SCRIPT trainval --experiment-dir VAE_SBF2 --exp-name VAE #--config $BASE/miann_ana/params/experiment_params_SBF2.py
#python $SCRIPT evaluate --experiment-dir VAE_all #VAE_SBF2
