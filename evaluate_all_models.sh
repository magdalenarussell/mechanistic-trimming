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
PRODUCTIVITY=$8

for LEFT_MOTIF_COUNT in {1..6..1}; do
    RIGHT_MOTIF_COUNT=0
    for MODEL_TYPE in distance_terminal_melting terminal_melting terminal_gc_content distance_terminal_gc_content; do
        COMMAND="sbatch -t 120 -c $NCPU -p $PARTITION -q $PARTITION $HOME/exomotif/evaluate_single_model.sh $ANNOTATION_TYPE $TRIM_TYPE $PRODUCTIVITY $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $MODEL_TYPE"
        $COMMAND
        echo "Running \`$COMMAND\`"
    done
    for RIGHT_MOTIF_COUNT in {1..6..1}; do
        for MODEL_TYPE in motif motif_distance motif_distance_terminal_gc_content motif_distance_terminal_melting motif_terminal_melting motif_terminal_gc_content; do
            COMMAND="sbatch -t 120 -c $NCPU -p $PARTITION -q $PARTITION $HOME/exomotif/evaluate_single_model.sh $ANNOTATION_TYPE $TRIM_TYPE $PRODUCTIVITY $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $MODEL_TYPE"
            $COMMAND
            echo "Running \`$COMMAND\`"
        done
        for LEFT_SIDE_TERMINAL_MELT_LENGTH in 5 10 15 20; do
            for MODEL_TYPE in motif_two_side_terminal_melting motif_distance_two_side_terminal_melting; do
                COMMAND="sbatch -t 120 -c $NCPU -p $PARTITION -q $PARTITION $HOME/exomotif/evaluate_single_model.sh $ANNOTATION_TYPE $TRIM_TYPE $PRODUCTIVITY $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $MODEL_TYPE $LEFT_SIDE_TERMINAL_MELT_LENGTH"
                $COMMAND
                echo "Running \`$COMMAND\`"
            done
        done
    done
done

RIGHT_MOTIF_COUNT=0
LEFT_MOTIF_COUNT=0
for MODEL_TYPE in distance; do
    COMMAND="sbatch -t 120 -c $NCPU -p $PARTITION -q $PARTITION $HOME/exomotif/evaluate_single_model.sh $ANNOTATION_TYPE $TRIM_TYPE $PRODUCTIVITY $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $MODEL_TYPE $LEFT_SIDE_TERMINAL_MELT_LENGTH"
    $COMMAND
    echo "Running \`$COMMAND\`"
done

for MODEL_TYPE in two_side_terminal_melting distance_two_side_terminal_melting; do
    for LEFT_SIDE_TERMINAL_MELT_LENGTH in 5 10 15 20; do
        COMMAND="sbatch -t 120 -c $NCPU -p $PARTITION -q $PARTITION $HOME/exomotif/evaluate_single_model.sh $ANNOTATION_TYPE $TRIM_TYPE $PRODUCTIVITY $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $MODEL_TYPE $LEFT_SIDE_TERMINAL_MELT_LENGTH"
        $COMMAND
        echo "Running \`$COMMAND\`"
    done
done

