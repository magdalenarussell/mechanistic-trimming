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

ANNOTATION_TYPE <<- args[1]
stopifnot(ANNOTATION_TYPE %in% c('igor', 'parsimony'))

TRIM_TYPE <<- args[2]
trim_types = list.files(path = 'scripts/gene_specific_functions/')
trim_types = str_sub(trim_types, end = -3)
stopifnot(TRIM_TYPE %in% trim_types)

MOTIF_TYPE <<- args[3] 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[4])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[5])
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[6])

UPPER_TRIM_BOUND <<- as.numeric(args[7]) 

TCR_REPERTOIRE_DATA <<- get(paste0('TCR_REPERTOIRE_DATA_', ANNOTATION_TYPE))

source('scripts/data_compilation_functions.R')

compile_all_motifs(TCR_REPERTOIRE_DATA)


