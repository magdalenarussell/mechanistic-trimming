#!/bin/bash

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate exomotif2 
set -eu

ANNOTATION_TYPE=$1
TRIM_TYPE=$2
MOTIF_TYPE=$3
GENE_WEIGHT_TYPE=$4
LEFT_MOTIF_COUNT=$5
RIGHT_MOTIF_COUNT=$6
UPPER_TRIM_BOUND=$7
PARTITION=$8
NCPU=$9
LEFT_SIDE_TERMINAL_MELT_LENGTH=${10}

for MODEL_TYPE in motif motif_distance distance; do
    COMMAND="sbatch -t 120 -c $NCPU -p $PARTITION -q $PARTITION $HOME/exomotif/evaluate_single_model.sh $ANNOTATION_TYPE $TRIM_TYPE $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $MODEL_TYPE"
    $COMMAND
    echo "Running \`$COMMAND\`"
done

for MODEL_TYPE in motif_two_side_terminal_melting motif_distance_two_side_terminal_melting two_side_terminal_melting distance_two_side_terminal_melting; do
    COMMAND="sbatch -t 120 -c $NCPU -p $PARTITION -q $PARTITION $HOME/exomotif/evaluate_single_model.sh $ANNOTATION_TYPE $TRIM_TYPE $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $MODEL_TYPE $LEFT_SIDE_TERMINAL_MELT_LENGTH"
    $COMMAND
    echo "Running \`$COMMAND\`"
done
