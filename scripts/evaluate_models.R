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

# NOTE: This method is only applicable for models fit across all subjects!
MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[7]
weight_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/sampling_procedure_functions/'))
weight_types = str_sub(weight_types, end = -3)
stopifnot(GENE_WEIGHT_TYPE %in% weight_types)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[8])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[9])

UPPER_TRIM_BOUND <<- as.numeric(args[10]) 
LOWER_TRIM_BOUND <<- as.numeric(args[11])

MODEL_TYPE <<- args[12]

TYPE <<- args[13]

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE) | grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[14])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

LOSS_GENE_WEIGHT <<- 'p_gene_given_subject'

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_evaluation_functions.R'))

# Compile data for all subjects
motif_data = aggregate_all_subject_data()

# Write model fit across all data
# fit_model_by_group(motif_data)

# compute loss
loss = evaluate_loss(motif_data)

# Append results to the end of a file
write_result_dt(loss$loss, rep(TYPE, length(loss$loss)), loss$model_parameter_count, loss$held_out_cluster_number, loss$held_out_genes)

if (TYPE == 'expected_log_loss'){
    file_name = get_per_run_model_evaluation_file_name(TYPE)
    file_name = str_replace(file_name, '/expected_log_loss/', '/expected_log_loss/raw/')

    losses = compile_result(loss$vect, TYPE, loss$model_parameter_count, loss$held_out_genes, held_out_clusters = NA, validation_gene_weighting = NA)

    fwrite(losses, file_name, sep = '\t')
}
