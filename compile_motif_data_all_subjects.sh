#!/bin/bash
#SBATCH --mail-user=magruss@uw.edu
#SBATCH --mail-type=END,FAIL
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate exomotif 
set -eu

Rscript $PWD/scripts/compile_motif_data.R $1 $2
