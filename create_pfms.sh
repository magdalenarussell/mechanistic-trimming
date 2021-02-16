#!/bin/bash
#SBATCH --mail-user=magruss@uw.edu
#SBATCH --mail-type=END,FAIL
source /home/mrussel2/miniconda3/etc/profile.d/conda.sh
conda activate exomotif 
set -eu

Rscript scripts/create_pfms.R $1 $2 $3 $4
