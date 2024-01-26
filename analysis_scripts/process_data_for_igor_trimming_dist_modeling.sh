#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming_jax
set -eu

MOD_OUTPUT_PATH=$1

param_path=$(python $PWD/mechanistic-trimming/analysis_scripts/process_for_igor_ligation-mh_model/get_igor_trimming_dist_params.py $MOD_OUTPUT_PATH)
echo "finished getting baseline igor parameters"

Rscript $PWD/mechanistic-trimming/analysis_scripts/process_for_igor_ligation-mh_model/convert_igor_trimming_dist_params.py $param_path
echo "finish preparing datasets"

