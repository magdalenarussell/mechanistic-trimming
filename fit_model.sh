#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming 
set -eu

ANNOTATION_TYPE=$1
DATA_GROUP=$2
TRIM_TYPE=$3
PRODUCTIVITY=$4
MOTIF_TYPE=$5
NCPU=$6
MODEL_GROUP=$7
GENE_WEIGHT_TYPE=$8
LEFT_MOTIF_COUNT=$9
RIGHT_MOTIF_COUNT=${10}
UPPER_TRIM_BOUND=${11}
LOWER_TRIM_BOUND=${12}
MODEL_TYPE=${13}
LEFT_SIDE_TERMINAL_MELT_LENGTH=${14:-NA}

Rscript $PWD/mechanistic-trimming/scripts/fit_model.R $ANNOTATION_TYPE $DATA_GROUP $TRIM_TYPE $PRODUCTIVITY $MOTIF_TYPE $NCPU $MODEL_GROUP $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $LOWER_TRIM_BOUND $MODEL_TYPE $LEFT_SIDE_TERMINAL_MELT_LENGTH 

