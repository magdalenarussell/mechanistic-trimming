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

GENE_WEIGHT_TYPE <<- 'p_gene_given_subject'

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

#TODO change this variable to a "group" variable and have data aggregation functions and model specific functions there...
REGRESSION_TYPE <<- 'all_subject'
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 4
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 4

UPPER_TRIM_BOUND <<- 18
LOWER_TRIM_BOUND <<- RIGHT_NUC_MOTIF_COUNT - 2 

source('scripts/data_compilation_functions.R')

#TODO add ability to look at different groups or individuals
motif_data = aggregate_all_subject_data()

#TODO next, model fitting

#TODO then, save coeffiecients, save model object, predict trimming curves...
