#!/bin/bash

ANNOTATION_TYPE=$1
TRIM_TYPE=$2
MOTIF_TYPE=$3
NCPU=$4
MODEL_GROUP=$5
GENE_WEIGHT_TYPE=$6
UPPER_TRIM_BOUND=$7

for MODEL in "motif" "motif_distance" "motif_distance_terminal_melting" "motif_terminal_melting" 
do 
    for LEFT in {2..6}
    do 
        for RIGHT in {2..6}
        do
            COMMAND="sbatch -c $NCPU fit_model.sh $ANNOTATION_TYPE $TRIM_TYPE $MOTIF_TYPE $NCPU $MODEL_GROUP $GENE_WEIGHT_TYPE $LEFT $RIGHT $UPPER_TRIM_BOUND $MODEL"
            echo $COMMAND
            $COMMAND
        done
    done
done

for MODEL in "distance" "terminal_melting" "distance_terminal_melting"
do
    COMMAND="sbatch -c $NCPU fit_model.sh $ANNOTATION_TYPE $TRIM_TYPE $MOTIF_TYPE $NCPU $MODEL_GROUP $GENE_WEIGHT_TYPE 0 0 $UPPER_TRIM_BOUND $MODEL"
    echo $COMMAND
    $COMMAND
done

