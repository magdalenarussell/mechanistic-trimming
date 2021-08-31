#!/bin/bash
#SBATCH --mail-user=magruss@uw.edu
#SBATCH --mail-type=END,FAIL

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate exomotif2 
set -eu

TRIM_TYPE=$1
MOTIF_TYPE=$2
NCPU=$3
GENE_WEIGHT_TYPE=$4
LEFT_MOTIF_COUNT=$5
RIGHT_MOTIF_COUNT=$6
UPPER_TRIM_BOUND=$7
MODEL_TYPE=$8

Rscript $PWD/scripts/evaluate_models.R $TRIM_TYPE $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $MODEL_TYPE 
