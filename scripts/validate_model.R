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

ANNOTATION_TYPE <<- args[1]
stopifnot(ANNOTATION_TYPE %in% c('igor', 'parsimony'))

TRIM_TYPE <<- args[2]
trim_types = list.files(path = 'scripts/gene_specific_functions/')
trim_types = str_sub(trim_types, end = -3)
stopifnot(TRIM_TYPE %in% trim_types)

PRODUCTIVITY <<- args[3]

MOTIF_TYPE <<- args[4] 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[5])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- args[6]

GENE_WEIGHT_TYPE <<- args[7]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[8])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[9])

UPPER_TRIM_BOUND <<- as.numeric(args[10]) 
LOWER_TRIM_BOUND <<- 2

MODEL_TYPE <<- args[11]

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[12])

if (!(grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE))){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

VALIDATION_DATA_DIR <<- args[13]
VALIDATION_TYPE <<- args[14]
VALIDATION_TRIM_TYPE <<- args[15]
VALIDATION_PRODUCTIVITY <<- args[16]
VALIDATION_GENE_NAME <<- paste0(substring(VALIDATION_TRIM_TYPE, 1, 1), '_gene')
stopifnot(VALIDATION_TYPE %in% c('validation_data_alpha', 'validation_data_beta', 'validation_data_gamma', 'validation_data_igh', 'validation_data_igk', 'validation_data_igl'))

source('scripts/model_fitting_functions.R')

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

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('scripts/model_evaluation_functions.R')

validation_data = aggregate_validation_data(directory = VALIDATION_DATA_DIR)

# compute loss
loss = evaluate_loss(validation_data, model)

# Append results to the end of a file
write_result_dt(loss$loss, rep(VALIDATION_TYPE, length(loss$loss)), loss$model_parameter_count, loss$held_out_cluster_number, loss$held_out_genes)
