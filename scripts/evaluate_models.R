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

# NOTE: This method is only applicable for models fit across all subjects!
MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[4]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[5])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[6])

UPPER_TRIM_BOUND <<- as.numeric(args[7]) 

MODEL_TYPE <<- args[8]

RESIDUAL_COMPARE_FEATURE <<- NULL

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('scripts/model_evaluation_functions.R')
source('plotting_scripts/residual_comparison_functions.R')

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

# Write model fit across all data
fit_model_by_group(motif_data)

predicted_trims = get_predicted_distribution_data()
per_gene_per_trim_resid = evaluate_model_per_gene(predicted_trims, type = 'per_gene_per_trim')
per_gene_resid = evaluate_model_per_gene(predicted_trims, type = 'per_gene')

# Append results to the end of a file
write_result_dt(log_loss, 'log_loss')
write_result_dt(per_gene_per_trim_resid, 'per_gene_per_trim')
write_result_dt(per_gene_resid, 'per_gene')
