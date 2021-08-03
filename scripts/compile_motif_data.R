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

NCPU <<- as.numeric(args[2])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[3])
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[4])

UPPER_TRIM_BOUND <<- 18
LOWER_TRIM_BOUND <<- 2
# TODO: Right now, all motifs have a bound of trim=2! This needs to change, but
# it is good to have consistent trim as we are comparing models
# LOWER_TRIM_BOUND <<- max(RIGHT_NUC_MOTIF_COUNT - 2, 1) 

source('scripts/data_compilation_functions.R')
print(LEFT_NUC_MOTIF_COUNT)
print(RIGHT_NUC_MOTIF_COUNT)
compile_all_motifs(TCR_REPERTOIRE_DATA)


