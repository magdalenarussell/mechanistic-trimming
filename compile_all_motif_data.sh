#!/bin/bash
set -eu

ANNOTATION_TYPE=$1
TRIM_TYPE=$2
MOTIF_TYPE=$3
NCPU=$4
UPPER_TRIM_BOUND=$5

for RIGHT in {0..6}; do
    for LEFT in {1..6}; do
        COMMAND="sbatch -c $NCPU compile_motif_data.sh $ANNOTATION_TYPE $TRIM_TYPE $MOTIF_TYPE $NCPU $LEFT $RIGHT $UPPER_TRIM_BOUND"
        echo $COMMAND
        $COMMAND
    done
done

COMMAND="sbatch -c $NCPU compile_motif_data.sh $ANNOTATION_TYPE $TRIM_TYPE $MOTIF_TYPE $NCPU 0 0 $UPPER_TRIM_BOUND"
echo $COMMAND
$COMMAND