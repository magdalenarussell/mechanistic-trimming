#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming 
set -eu

ANNOTATION_TYPE=$1
PARAM_GROUP=$2
NCPU=$3
GENE_WEIGHT_TYPE=$4
LEFT_MOTIF_COUNT=$5
RIGHT_MOTIF_COUNT=$6
MODEL_TYPE=$7

Rscript $PWD/mechanistic-trimming/scripts/fit_model.R $ANNOTATION_TYPE $PARAM_GROUP $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $MODEL_TYPE

