source('config/config.R')

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

TRIM_TYPE <<- args[1]
stopifnot(TRIM_TYPE == 'v_trim')

MOTIF_TYPE <<- args[2] 
stopifnot(MOTIF_TYPE %in% c('bounded', 'unbounded', 'unbounded_no_pnuc'))

NCPU <<- as.numeric(args[3])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- args[4]

GENE_WEIGHT_TYPE <<- args[5]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[6])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[7])

UPPER_TRIM_BOUND <<- as.numeric(args[8]) 

MODEL_TYPE <<- args[9]

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')

# Compile data for all subjects
motif_data = aggregate_all_subject_data()

# Fit model, write predicted distribution and pwm files
fit_model_by_group(motif_data)
