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
stopifnot(MOTIF_TYPE %in% c('bounded', 'unbounded'))

NCPU <<- as.numeric(args[3])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[4]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[5])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[6])

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('scripts/model_evaluation_functions.R')

# Compile data for all subjects
motif_data = aggregate_all_subject_data()

# Generate a held out sample and motif data subset
sample_data = generate_hold_out_sample(motif_data, sample_size = length(unique(motif_data$gene)))
motif_data_subset = sample_data$motif_data_subset
sample = sample_data$sample

# Fit model to the motif_data_subset
model = fit_model(motif_data_subset)

# Compute conditional logistic loss value for held out sample using model
log_loss = calculate_cond_log_loss(model, sample)

# Append results to the end of a file
write_result_dt(log_loss, sample)
    
