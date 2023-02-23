#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming 
set -eu

ANNOTATION_TYPE=$1
TRIM_TYPE=$2
PRODUCTIVITY=$3
MOTIF_TYPE=$4
NCPU=$5
GENE_WEIGHT_TYPE=$6
LEFT_MOTIF_COUNT=$7
RIGHT_MOTIF_COUNT=$8
UPPER_TRIM_BOUND=$9
MODEL_TYPE=${10}
PROP=${11}
LEFT_SIDE_TERMINAL_MELT_LENGTH=${12:-NA}

Rscript $PWD/scripts/subsampling_experiment.R $ANNOTATION_TYPE $TRIM_TYPE $PRODUCTIVITY $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $MODEL_TYPE $PROP $LEFT_SIDE_TERMINAL_MELT_LENGTH 

