#!/bin/bash
#SBATCH --mail-user=magruss@uw.edu
#SBATCH --mail-type=END,FAIL
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate exomotif2 
set -eu

ANNOTATION_TYPE=$1
TRIM_TYPE=$2
PRODUCTIVITY=$3
MOTIF_TYPE=$4
NCPU=$5
LEFT_MOTIF_COUNT=$6
RIGHT_MOTIF_COUNT=$7
UPPER_TRIM_BOUND=$8

Rscript $PWD/scripts/compile_motif_data.R $ANNOTATION_TYPE $TRIM_TYPE $PRODUCTIVITY $MOTIF_TYPE $NCPU $LEFT_MOTIF_COUNT $RIGHT_MOTIF_COUNT $UPPER_TRIM_BOUND
