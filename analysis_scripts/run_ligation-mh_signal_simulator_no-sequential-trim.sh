#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming_jax
set -eu

MOD_OUTPUT_PATH=$1
NCPU=$2
TRIM_SAMPLING_TYPE=$3
LIGATION_MH_PARAM=$4

param_path=$(python $PWD/mechanistic-trimming/analysis_scripts/ligation-mh_signal_simulator_scripts/get_igor_params.py $MOD_OUTPUT_PATH)
echo "finished getting baseline igor parameters"

Rscript $PWD/mechanistic-trimming/analysis_scripts/ligation-mh_signal_simulator_no-sequential-trim/process_data_for_ligation-mh_signal_simulation.R $NCPU $TRIM_SAMPLING_TYPE $LIGATION_MH_PARAM
