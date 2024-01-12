#!/bin/bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mechanistic-trimming_jax
set -eu

NCPU=$1

path=$(Rscript $PWD/mechanistic-trimming/analysis_scripts/process_data_for_negative_ligation-mh_signal_simulation.R)
echo "finished processing data for model prediction"

python $PWD/mechanistic-trimming/jax_scripts/predict.py $path igor_sim_alpha nonproductive_v-j_trim_ligation-mh 1 2 motif_two-side-base-count-beyond False $NCPU pre_simulation
echo "finished making predictions"

Rscript $PWD/mechanistic-trimming/analysis_scripts/simulate_data_for_negative_ligation-mh_signal_simulation.R
echo "finished simulating data"

