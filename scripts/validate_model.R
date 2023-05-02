source('mechanistic-trimming/config/config.R')

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

ANNOTATION_TYPE <<- args[1]
stopifnot(ANNOTATION_TYPE %in% c('igor', 'parsimony', 'alpha'))

DATA_GROUP <<- args[2]
group_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/data_grouping_functions/'))
group_types = str_sub(group_types, end = -3)
stopifnot(DATA_GROUP %in% group_types)

TRIM_TYPE <<- args[3]
trim_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/gene_specific_functions/'))
trim_types = str_sub(trim_types, end = -3)
stopifnot(TRIM_TYPE %in% trim_types)

PRODUCTIVITY <<- args[4]

MOTIF_TYPE <<- args[5] 
motif_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/motif_class_functions/'))
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[6])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- args[7]

GENE_WEIGHT_TYPE <<- args[8]
weight_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/sampling_procedure_functions/'))
weight_types = str_sub(weight_types, end = -3)
stopifnot(GENE_WEIGHT_TYPE %in% weight_types)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[9])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[10])

UPPER_TRIM_BOUND <<- as.numeric(args[11]) 
LOWER_TRIM_BOUND <<- as.numeric(args[12])

MODEL_TYPE <<- args[13]

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[14])

if (!(grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE))){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

VALIDATION_DATA_DIR <<- args[15]
VALIDATION_TYPE <<- args[16]
VALIDATION_TRIM_TYPE <<- args[17]
VALIDATION_PRODUCTIVITY <<- args[18]
VALIDATION_GENE_NAME <<- paste0(substring(VALIDATION_TRIM_TYPE, 1, 1), '_gene')
stopifnot(VALIDATION_TYPE %in% c('validation_data_alpha', 'validation_data_beta', 'validation_data_gamma', 'validation_data_delta', 'validation_data_igh', 'validation_data_igk', 'validation_data_igl', 'igor'))

LOSS_GENE_WEIGHT <<- args[19]
stopifnot(LOSS_GENE_WEIGHT %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform', 'p_gene_marginal_all_seqs', 'p_gene_given_subject_all_seqs'))

source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))

if (MODEL_TYPE != 'null') {
    model = load_model()
} else {
    model = 'null'
} 

ANNOTATION_TYPE <<- VALIDATION_TYPE
TYPE <<- 'validation_data' 
TRIM_TYPE <<- VALIDATION_TRIM_TYPE
GENE_NAME <<- VALIDATION_GENE_NAME
PRODUCTIVITY <<- VALIDATION_PRODUCTIVITY

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_evaluation_functions.R'))

validation_data = aggregate_validation_data(directory = VALIDATION_DATA_DIR)

# compute loss
loss = evaluate_loss(validation_data, model)

# Append results to the end of a file
write_result_dt(loss$loss, rep(VALIDATION_TYPE, length(loss$loss)), loss$model_parameter_count, loss$held_out_cluster_number, loss$held_out_genes, validation_gene_weighting = LOSS_GENE_WEIGHT)
