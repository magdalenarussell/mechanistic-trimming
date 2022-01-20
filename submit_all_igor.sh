#!/bin/bash

set -eu

RAW_FILE_PATH=$1
TEMP_DIR=$2
OUTPUT_DIR=$3
ANNOTATION_COUNT=$4
NCPU=$5

for RAW_FILE in $RAW_FILE_PATH/*.tsv; do
    COMMAND="sbatch -c $NCPU igor_scripts/run_igor_all.sh $RAW_FILE $TEMP_DIR $OUTPUT_DIR $ANNOTATION_COUNT $NCPU"
    echo $COMMAND
    $COMMAND
done  

