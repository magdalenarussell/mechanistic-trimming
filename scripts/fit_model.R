source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)


args = commandArgs(trailingOnly=TRUE)

TRIM_TYPE <<- args[1]
stopifnot(TRIM_TYPE == 'v_trim')

MOTIF_TYPE <<- 'bounded'
stopifnot(MOTIF_TYPE %in% c('bounded', 'unbounded'))

NCPU <<- args[2]

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- args[3]

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 4
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 4

UPPER_TRIM_BOUND <<- 18
LOWER_TRIM_BOUND <<- RIGHT_NUC_MOTIF_COUNT - 2 

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')

# Compile data for all subjects
motif_data = aggregate_all_subject_data()

# TODO add functionality to compare between individuals, groups

#TODO next, model fitting

#TODO then, save coeffiecients, save model object, predict trimming curves...
