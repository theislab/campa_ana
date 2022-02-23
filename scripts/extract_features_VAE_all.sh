#! /bin/bash

#SBATCH -o features_%j.out
#SBATCH -e features_%j.out
#SBATCH -J features
#SBATCH -p cpu_p
#SBATCH -c 8
#SBATCH --mem 127G
#SBATCH -t 40:00:00
#SBATCH --nice=10000

source ~/.bashrc
conda activate pelkmans-3.9

BASE=/home/icb/hannah.spitzer/projects/pelkmans/software_new
SCRIPT=$BASE/campa/cli/extract_features.py


# -- VAE_all experiment --
# intensity
#python $SCRIPT VAE_all/CondVAE_pert-CC "clustering_res0.5" --save-name features.h5ad intensity
#python $SCRIPT VAE_all/CondVAE_pert-CC "clustering_res0.5" --cluster-col annotation --save-name features_annotation.h5ad intensity

# co-occurrence (for annotated data)
#python $SCRIPT \
#    VAE_all/CondVAE_pert-CC "clustering_res0.5" intensity co-occurrence \
#    --cluster-col annotation --save-name features_annotation.h5ad --co-logspace --co-nsteps 20 --num-processes 8 \
#    --data-dir \
#    "184A1_unperturbed/I09" "184A1_unperturbed/I11" "184A1_unperturbed/J10" "184A1_unperturbed/J12" \
#    "184A1_DMSO/I14" "184A1_DMSO/J16" "184A1_AZD4573/I13" "184A1_AZD4573/I17" "184A1_AZD4573/J14" \
#    "184A1_AZD4573/J18" "184A1_AZD4573/J21" "184A1_CX5461/I18" "184A1_CX5461/J09" "184A1_CX5461/J22" \
#    "184A1_TSA/I16" "184A1_TSA/J13" "184A1_TSA/J20" "184A1_triptolide/I10" "184A1_triptolide/J15" \
#    "184A1_meayamycin/I12" "184A1_meayamycin/I20" 

# object stats (for annotated data)
python $SCRIPT VAE_all/CondVAE_pert-CC "clustering_res0.5" --cluster-col annotation --save-name features_annotation.h5ad object-stats \
    --stats-features 'area' 'circularity' 'elongation' 'extent'

    