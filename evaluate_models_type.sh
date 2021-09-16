#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate exomotif2 
set -eu

TRIM_TYPE=$1
MOTIF_TYPE=$2
GENE_WEIGHT_TYPE=$3
UPPER_TRIM_BOUND=$4
PARTITION=$5
NCPU=$6
MODEL_TYPE=$7

for LEFT_MOTIF_COUNT in {2..6..1}; do
    for RIGHT_MOTIF_COUNT in {2..6..1}; do
        COMMAND="sbatch -c $NCPU -p $PARTITION -q $PARTITION $HOME/exomotif/evaluate_single_model.sh $TRIM_TYPE $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $MODEL_TYPE"
        $COMMAND
        echo "Running \`$COMMAND\`"
    done
done
