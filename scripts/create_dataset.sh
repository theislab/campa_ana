#! /bin/bash

#SBATCH -o create_dataset_%j.out
#SBATCH -e create_dataset_%j.out
#SBATCH -J create_dataset
#SBATCH -p cpu_p
#SBATCH -c 1
#SBATCH --mem-per-cpu 89G
#SBATCH -t 8:00:00
#SBATCH --nice=10000

# set up
#source ~/.bashrc
#conda activate pelkmans-3.9

# run script
CUR_DIR=$(dirname $0)
SCRIPT=$CUR_DIR/../cli/create_dataset.py

python $SCRIPT $CUR_DIR/params/data_params_all.py

