library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)

TRIM_TYPE <<- 'v_trim' 
stopifnot(TRIM_TYPE == 'v_trim')

MOTIF_TYPE <<- 'pnuc_motif' 
stopifnot(MOTIF_TYPE %in% c('pnuc_motif', 'no_pnuc_motif'))

PFM_TYPE <<- 'unbounded'
stopifnot(PFM_TYPE %in% c('unbounded'))

NCPU <<- 1

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 4
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 4

source('scripts/pfm_functions.R')

PFM_no_pnucs = as.matrix(fread('_ignore/pfms/complete/ALL_unbounded_no_pnuc_motif_ALL.tsv'), rownames = 1)
PWM_no_pnucs = calculate_PWM_from_PFM(PFM_no_pnucs)
