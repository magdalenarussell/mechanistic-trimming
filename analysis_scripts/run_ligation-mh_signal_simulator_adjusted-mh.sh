#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming_jax
set -eu

MOD_OUTPUT_PATH=$1
PARAM_GROUP=$2
NCPU=$3
TRIM_SAMPLING_TYPE=$4
LIGATION_MH_PARAM=$5

param_path=$(python $PWD/mechanistic-trimming/analysis_scripts/ligation-mh_signal_simulator/get_igor_params.py $MOD_OUTPUT_PATH)
echo "finished getting baseline igor parameters"

Rscript $PWD/mechanistic-trimming/analysis_scripts/ligation-mh_signal_simulator_adjusted-mh/process_data_for_ligation-mh_signal_simulation.R $PARAM_GROUP $NCPU $TRIM_SAMPLING_TYPE $LIGATION_MH_PARAM