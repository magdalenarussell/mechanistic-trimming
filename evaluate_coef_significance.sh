#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming 
set -eu

ANNOTATION_TYPE=$1
DATA_GROUP=$2
PARAM_GROUP=$3
NCPU=$4
GENE_WEIGHT_TYPE=$5
LEFT_MOTIF_COUNT=$6
RIGHT_MOTIF_COUNT=$7
MODEL_TYPE=$8

Rscript $PWD/mechanistic-trimming/scripts/evaluate_coeffs.R $ANNOTATION_TYPE $DATA_GROUP $PARAM_GROUP $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $MODEL_TYPE 

