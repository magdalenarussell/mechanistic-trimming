TRIM_TYPE <<- 'v_trim' 
stopifnot(TRIM_TYPE == 'v_trim')

MOTIF_TYPE <<- 'pnuc_motif' 

PFM_TYPE <<- 'bounded'

NCPU <<- 4 
REGRESSION_TYPE <<- 'all_subject' 
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 4
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 4
UPPER_TRIM_BOUND <<- 18
LOWER_TRIM_BOUND <<- 2


library(speedglm)
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(NCPU)
library(Biostrings)
library(RhpcBLASctl)
omp_set_num_threads(NCPU)
blas_set_num_threads(NCPU)

library(speedglm)
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(NCPU)
library(Biostrings)
library(RhpcBLASctl)
omp_set_num_threads(NCPU)
blas_set_num_threads(NCPU)

source('scripts/pfm_functions.R')

OUTPUT_PATH = '/fh/fast/matsen_e/shared/tcr-gwas/exomotif/motif_data'

together = concatenate_motifs_all_subjects()


