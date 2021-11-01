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
#source ~/.bashrc
#module load Anaconda3

# run script
CUR_DIR=$(dirname $0)
SCRIPT=$CUR_DIR/../cli/download_mpp.py
PARAMS=$CUR_DIR/params/download_params.py

echo "running $SCRIPT from $CUR_DIR"
python $SCRIPT --params $PARAMS --key unperturbed
#python $SCRIPT --params $PARAMS --key AZD4573 --experiment_name $EXPERIMENT
# TODO add remaining calls
#python $BASEDIR/../cli/download_mpp.py -w I13 I17 J14 J18 J21 -o ~/NascentRNA/local_data/184A1_hannah_AZD4573
#python download_mpp.py ~/NascentRNA/tm_credentials -w I18 J09 J22 -o ~/NascentRNA/local_data/184A1_hannah_CX5461
#python download_mpp.py ~/NascentRNA/tm_credentials -w I14 J16 -o ~/NascentRNA/local_data/184A1_hannah_DMSO
#python download_mpp.py ~/NascentRNA/tm_credentials -w I16 J13 J20 -o ~/NascentRNA/local_data/184A1_hannah_TSA
#python download_mpp.py ~/NascentRNA/tm_credentials -w I10 J15 -o ~/NascentRNA/local_data/184A1_hannah_triptolide
#python download_mpp.py ~/NascentRNA/tm_credentials -w H07 I08 J07 -o ~/NascentRNA/local_data/184A1_hannah_EU10
#python download_mpp.py ~/NascentRNA/tm_credentials -w I12 I20 -o ~/NascentRNA/local_data/184A1_hannah_meayamycin

#python download_mpp.py ~/NascentRNA/tm_credentials -w K18 L19 M18 -o ~/NascentRNA/local_data/184A1_hannah_SBF2 --object Cells
#python download_mpp.py ~/NascentRNA/tm_credentials -w K19 L18 M19 -o ~/NascentRNA/local_data/184A1_hannah_scrambled --object Cells