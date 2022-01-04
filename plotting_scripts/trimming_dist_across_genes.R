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
trim_types = list.files(path = 'scripts/gene_specific_functions/')
trim_types = str_sub(trim_types, end = -3)
stopifnot(TRIM_TYPE %in% trim_types)

PRODUCTIVITY <<- args[2]

NCPU <<- as.numeric(args[3])

MOTIF_TYPE <<- 'bounded'
MODEL_GROUP <<- 'all_subjects'

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

source('scripts/data_compilation_functions.R')
source('plotting_scripts/plotting_functions.R')
compile_all_motifs(TCR_REPERTOIRE_DATA)


