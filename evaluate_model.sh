#!/bin/bash
#SBATCH --mail-user=magruss@uw.edu
#SBATCH --mail-type=END,FAIL
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate exomotif2 
set -eu

Rscript $PWD/scripts/evaluate_models.R $1 $2 $3 $4 $5 $6 $7 
