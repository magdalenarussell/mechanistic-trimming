source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(mclogit)
library(matrixcalc)
library(rdist)
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

MODEL_GROUP <<- args[3]

GENE_WEIGHT_TYPE <<- args[4]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[5])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[6])

GENE_PROB_TYPE <<- args[7]
stopifnot(GENE_PROB_TYPE %in% c('uniform', 'non_uniform'))

UPPER_TRIM_BOUND <<- 18
LOWER_TRIM_BOUND <<- RIGHT_NUC_MOTIF_COUNT - 2 

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('scripts/simulation_functions.R')

SIM_OUTPUT_LOC <<- '/fh/scratch/delete30/matsen_e/mrussel2/motif_sims'
PWM = generate_random_PWM(gene_count = 10) 

gene_count = 10
subject_names = c('A', 'B', 'C', 'D', 'E')

# simulate data using given probability distribution for genes and PWM
sim_data = get_complete_sim_data(subject_names, gene_count, PWM, gene_probs = GENE_PROB_TYPE, gene_probs_type = GENE_PROB_TYPE)
    
# fit model
file_path = get_file_path(gene_count, gene_probs_type = GENE_PROB_TYPE)
motif_data = aggregate_all_subject_data(directory = file_path)
fit_PWM = fit_model_by_group(motif_data, write_coeffs = FALSE)

# calculate euclidean distance of coefficient matrix versus original PWM
fit_PWM_mat = as.matrix(fit_PWM[, -c('base', 'model_group')])
rownames(fit_PWM_mat) = fit_PWM$base
dist = cdist(fit_PWM_mat, PWM)
print(sum(dist))

# NOTE: I have found that, as expected, when we using the p_gene_marginal gene
# weight type for model fitting, the coefficient matrix is closer to the
# original PWM compared to when we use the p_gene_given_subject weight type for
# model fitting. This relationship is more pronounced gene choice originates
# from a non-uniform distribution for each subject.
