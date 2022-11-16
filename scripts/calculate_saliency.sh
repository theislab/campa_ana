#! /bin/bash

#SBATCH -o calculate_saliency_%j.out
#SBATCH -e calculate_saliency_%j.out
#SBATCH -J jupyterlab
#SBATCH -p cpu_p
#SBATCH -c 1
##SBATCH --mem=63G
#SBATCH --mem=127G
#SBATCH -t 23:00:00
#SBATCH --nice=10000
#SBATCH --exclude=ibis216-010-[005,007,068-071]

# set up
source ~/.bashrc
conda activate campa

# run script for calculating saliencies
BASE=/home/icb/hannah.spitzer/projects/pelkmans/software_new/campa_ana

python $BASE/figure_notebooks/calculate_saliency.py

