#! /bin/bash

#SBATCH -o features_%j.out
#SBATCH -e features_%j.out
#SBATCH -J features
#SBATCH -p cpu_p
#SBATCH -c 8
#SBATCH --mem 127G
#SBATCH -t 40:00:00
#SBATCH --nice=10000

# set up
source ~/.bashrc
conda activate pelkmans-3.9

# run script
BASE=/home/icb/hannah.spitzer/projects/pelkmans/software_new
SCRIPT=$BASE/campa/cli/extract_features.py

# calculate features for VAE_all experiment according to feature_params_all.py
python $SCRIPT $BASE/campa_ana/params/feature_params_all.py
# calculate features for VAE_SBF2 experiment according to feature_params_SBF2.py
# NOTE: calculation of co-occurrence scores takes up to 40 hours per well for this experiment
#python $SCRIPT $BASE/campa_ana/params/feature_params_SBF2.py
