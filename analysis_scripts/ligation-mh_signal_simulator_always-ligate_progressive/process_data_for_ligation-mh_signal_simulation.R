source('mechanistic-trimming/config/config.R')
source('mechanistic-trimming/config/file_paths.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- as.numeric(args[1])
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_ligation-mh'

TRIMMING_PROB_TYPE <<- 'uniform'

LIGATION_MH_PARAM <<- 0

source(paste0(MOD_PROJECT_PATH,'/analysis_scripts/ligation-mh_signal_simulator_scripts/ligation-mh_simulator_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

old_ANNOTATION_TYPE <<- ANNOTATION_TYPE
ANNOTATION_TYPE <<- paste0(old_ANNOTATION_TYPE, '_from_', TRIMMING_PROB_TYPE, '_MHprob', LIGATION_MH_PARAM)

filename = processed_data_path()
data = fread(filename)

MODEL_TYPE <<- 'ligation-mh'
LEFT_NUC_MOTIF_COUNT <<- 0
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 0 

filename = processed_data_path()
fwrite(data, filename, sep = '\t')

# get gene pair frequencies
freq = data[, sum(count), by = v_gene][order(-V1)]

## run experiment to add gene pairs in to the model training set progressively
### start with the most frequently used gene pair
prog_train = data.table()
pairs = 0
for (i in seq(nrow(freq))){
    # add new pair to training dataset
    v = freq[i]$v_gene
    temp = data[v_gene == v]
    prog_train = rbind(prog_train, temp)

    v = str_replace(v, '\\/', '.')
    output_path = file.path(dirname(filename), paste0('progressively_added_data_', v, '_', pairs, '.tsv'))

    fwrite(prog_train, output_path, sep = '\t')
    pairs = pairs + 1
}
