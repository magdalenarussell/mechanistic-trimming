#!/bin/bash

ANNOTATION_TYPE=$1
TRIM_TYPE=$2
MOTIF_TYPE=$3
NCPU=$4
MODEL_GROUP=$5
GENE_WEIGHT_TYPE=$6
UPPER_TRIM_BOUND=$7

for LEFT in {2..6}
do 
    for RIGHT in {2..6}
    do
        for MODEL in "motif" "motif_distance" "motif_distance_terminal_melting" "motif_terminal_melting" "motif_distance_terminal_melting_NN" "motif_terminal_melting_NN" 
        do 
            COMMAND="sbatch -c $NCPU fit_model.sh $ANNOTATION_TYPE $TRIM_TYPE $MOTIF_TYPE $NCPU $MODEL_GROUP $GENE_WEIGHT_TYPE $LEFT $RIGHT $UPPER_TRIM_BOUND $MODEL"
            echo $COMMAND
            $COMMAND
        done
    done
done

for MODEL in "distance" "terminal_melting" "distance_terminal_melting" "terminal_melting_NN" "distance_terminal_melting_NN"
do
    COMMAND="sbatch -c $NCPU fit_model.sh $ANNOTATION_TYPE $TRIM_TYPE $MOTIF_TYPE $NCPU $MODEL_GROUP $GENE_WEIGHT_TYPE 0 0 $UPPER_TRIM_BOUND $MODEL"
    echo $COMMAND
    $COMMAND
done

for LEFT_SIDE_TERMINAL_MELT_LENGTH in 5 10 15 20; do 
    for MODEL in "two_side_terminal_melting" "distance_two_side_terminal_melting" "two_side_terminal_melting_NN" "distance_two_side_terminal_melting_NN"; do
        COMMAND="sbatch -c $NCPU fit_model.sh $ANNOTATION_TYPE $TRIM_TYPE $MOTIF_TYPE $NCPU $MODEL_GROUP $GENE_WEIGHT_TYPE 0 0 $UPPER_TRIM_BOUND $MODEL $LEFT_SIDE_TERMINAL_MELT_LENGTH"
        echo $COMMAND
        $COMMAND
    done
    for LEFT in {2..6}; do
        for RIGHT in {2..6}; do
            for MODEL in "motif_distance_two_side_terminal_melting" "motif_two_side_terminal_melting" "motif_distance_two_side_terminal_melting_NN" "motif_two_side_terminal_melting_NN"; do
                COMMAND="sbatch -c $NCPU fit_model.sh $ANNOTATION_TYPE $TRIM_TYPE $MOTIF_TYPE $NCPU $MODEL_GROUP $GENE_WEIGHT_TYPE $LEFT $RIGHT $UPPER_TRIM_BOUND $MODEL $LEFT_SIDE_TERMINAL_MELT_LENGTH"
                echo $COMMAND
                $COMMAND
            done
        done
    done
done
