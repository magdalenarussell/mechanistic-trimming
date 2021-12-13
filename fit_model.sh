#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate exomotif2 
set -eu

LEFT_SIDE_TERMINAL_MELT_LENGTH=${11:-NA}

Rscript $PWD/scripts/fit_model.R $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} $LEFT_SIDE_TERMINAL_MELT_LENGTH
