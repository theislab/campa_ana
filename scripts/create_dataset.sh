#! /bin/bash

#SBATCH -o create_dataset_%j.out
#SBATCH -e create_dataset_%j.out
#SBATCH -J create_dataset
#SBATCH -p cpu_p
#SBATCH -c 1
#SBATCH --mem-per-cpu 255G
#SBATCH -t 8:00:00
#SBATCH --nice=10000

# set up
source ~/.bashrc
conda activate campa

BASE=/home/icb/hannah.spitzer/projects/pelkmans/software_new
# run script
campa create_dataset $BASE/campa_ana/params/data_params_all.py
campa create_dataset $BASE/campa_ana/params/data_params_all_noneigh.py
campa create_dataset $BASE/campa_ana/params/data_params_SBF2.py

