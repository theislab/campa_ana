#! /bin/bash
#
#SBATCH --time=12:00:00
#SBATCH --mem=16G
#
### write both output and errors to file `run_mcu.NNN.log`
#SBATCH -o run_download_mpp.%j.log
#SBATCH -e run_download_mpp.%j.log
#

# set up
source ~/.bashrc
module load Anaconda3
conda activate pelkmans-3.9

# run script
CUR_DIR=$(dirname $0)
CUR_DIR=/data/homes/hannah/NascentRNA/software_new/campa_ana/scripts
SCRIPT=$CUR_DIR/download_mpp.py
#PARAMS=$CUR_DIR/../params/download_params.py
PARAMS=$CUR_DIR/../params/download_params_ilastik.py

echo "running $SCRIPT from $CUR_DIR"
python $SCRIPT --params $PARAMS --key unperturbed
python $SCRIPT --params $PARAMS --key AZD4573
python $SCRIPT --params $PARAMS --key CX5461
python $SCRIPT --params $PARAMS --key DMSO
python $SCRIPT --params $PARAMS --key TSA
python $SCRIPT --params $PARAMS --key triptolide
python $SCRIPT --params $PARAMS --key EU10
python $SCRIPT --params $PARAMS --key meayamycin
#python $SCRIPT --params $PARAMS --key SBF2
#python $SCRIPT --params $PARAMS --key scrambled