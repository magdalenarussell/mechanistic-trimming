#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate exomotif2 
set -eu

ANNOTATION_TYPE=$1
TRIM_TYPE=$2
MOTIF_TYPE=$3
GENE_WEIGHT_TYPE=$4
UPPER_TRIM_BOUND=$5
PARTITION=$6
NCPU=$7
MODEL_TYPE=$8

for LEFT_MOTIF_COUNT in {2..6..1}; do
    for RIGHT_MOTIF_COUNT in {2..6..1}; do
        COMMAND="sbatch -c $NCPU -p $PARTITION -q $PARTITION $HOME/exomotif/evaluate_single_model.sh $ANNOTATION_TYPE $TRIM_TYPE $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $MODEL_TYPE"
        $COMMAND
        echo "Running \`$COMMAND\`"
    done
done
