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
VALIDATION_DATA_DIR=$8
VALIDATION_TYPE=$9
VALIDATION_TRIM_TYPE=${10}
VALIDATION_PRODUCTIVITY=${11}

Rscript $PWD/mechanistic-trimming/scripts/analysis_scripts/evaluate_relative_importance.R $ANNOTATION_TYPE $DATA_GROUP $PARAM_GROUP $NCPU $MODEL_GROUP $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $VALIDATION_DATA_DIR $VALIDATION_TYPE $VALIDATION_TRIM_TYPE $VALIDATION_PRODUCTIVITY


