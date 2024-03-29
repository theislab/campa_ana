#! /bin/bash

#SBATCH -o start_notebook_%j.out
#SBATCH -e start_notebook_%j.out
#SBATCH -J jupyterlab
#SBATCH -p cpu_p
#SBATCH -c 1
##SBATCH --mem=63G
#SBATCH --mem=127G
#SBATCH -t 10:00:00
#SBATCH --nice=10000
#SBATCH --exclude=ibis216-010-[005,007,068-071]

# set up
source ~/.bashrc
conda activate campa

# run script
BASE=/home/icb/hannah.spitzer/projects/pelkmans/software_new
# prevent othher people from reading the log file
chmod og-r start_notebook_${SLURM_JOB_ID}.out
jupyter-lab --no-browser --ip 0.0.0.0 --notebook-dir $BASE
