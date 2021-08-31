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

NCPU <<- as.numeric(args[2])

MOTIF_TYPE <<- 'bounded'
MODEL_GROUP <<- 'all_subjects'

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

source('scripts/data_compilation_functions.R')
source('plotting_scripts/plotting_functions.R')
compile_all_motifs(TCR_REPERTOIRE_DATA)


