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

for LEFT_MOTIF_COUNT in {2..6..1}; do
    for MODEL_TYPE_1 in distance_terminal_gc_content distance_terminal_melting terminal_melting terminal_gc_content ; do
        COMMAND1="sbatch -t 120 -c $NCPU -p $PARTITION -q $PARTITION $HOME/exomotif/evaluate_single_model.sh $TRIM_TYPE $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT 0 $UPPER_TRIM_BOUND $MODEL_TYPE_1"
        $COMMAND1
        echo "Running \`$COMMAND1\`"
    done
    for RIGHT_MOTIF_COUNT in {2..6..1}; do
        for MODEL_TYPE in motif motif_distance motif_terminal_gc_content motif_distance_terminal_gc_content motif_terminal_melting motif_distance_terminal_melting ; do
            COMMAND="sbatch -t 120 -c $NCPU -p $PARTITION -q $PARTITION $HOME/exomotif/evaluate_single_model.sh $TRIM_TYPE $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND $MODEL_TYPE"
            $COMMAND
            echo "Running \`$COMMAND\`"
        done
    done
done

sbatch -t 120 -c $NCPU -p $PARTITION -q $PARTITION $HOME/exomotif/evaluate_single_model.sh $TRIM_TYPE $MOTIF_TYPE $NCPU $GENE_WEIGHT_TYPE 0 0 $UPPER_TRIM_BOUND distance
echo "Running distance only"
