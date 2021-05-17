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

MOTIF_TYPE <<- 'pnuc_motif' 
stopifnot(MOTIF_TYPE %in% c('pnuc_motif', 'no_pnuc_motif'))

PFM_TYPE <<- 'bounded'
stopifnot(PFM_TYPE %in% c('bounded', 'unbounded'))

NCPU <<- args[2]

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

REGRESSION_TYPE <<- 'all_subject'
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 4
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 4

UPPER_TRIM_BOUND <<- 18
LOWER_TRIM_BOUND <<- 2

source('scripts/pfm_functions.R')

compile_all_motifs(TCR_REPERTOIRE_DATA)

