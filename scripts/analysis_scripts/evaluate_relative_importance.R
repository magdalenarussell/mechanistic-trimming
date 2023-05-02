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
stopifnot(ANNOTATION_TYPE %in% c('igor', 'parsimony'))

DATA_GROUP <<- args[2]
group_types = list.files(path = 'scripts/data_grouping_functions/')
group_types = str_sub(group_types, end = -3)
stopifnot(DATA_GROUP %in% group_types)

TRIM_TYPE <<- args[3]
trim_types = list.files(path = 'mechanistic-trimming/scripts/gene_specific_functions/')
trim_types = str_sub(trim_types, end = -3)
stopifnot(TRIM_TYPE %in% trim_types)

PRODUCTIVITY <<- args[4]

MOTIF_TYPE <<- args[5] 
motif_types = list.files(path = 'mechanistic-trimming/scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[6])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- args[7]

GENE_WEIGHT_TYPE <<- args[8]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[9])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[10])

UPPER_TRIM_BOUND <<- as.numeric(args[11]) 
LOWER_TRIM_BOUND <<- as.numeric(args[12])

MODEL_TYPE <<- 'motif_two-side-base-count-beyond'

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[13])

if (!(grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE))){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

VALIDATION_DATA_DIR <<- args[14]
VALIDATION_TYPE <<- args[15]
VALIDATION_TRIM_TYPE <<- args[16]
VALIDATION_PRODUCTIVITY <<- args[17]
VALIDATION_GENE_NAME <<- paste0(substring(VALIDATION_TRIM_TYPE, 1, 1), '_gene')
stopifnot(VALIDATION_TYPE %in% c('validation_data_alpha', 'validation_data_beta', 'validation_data_gamma', 'validation_data_delta', 'validation_data_igh', 'validation_data_igk', 'validation_data_igl', 'igor'))

LOSS_GENE_WEIGHT <<- 'p_gene_given_subject' 

source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/analysis_scripts/rel_importance_functions.R'))

model = load_model()
pwm = get_model_coefficient_data()
output_file_name = get_per_run_rel_importance_model_file_name()

ANNOTATION_TYPE <<- VALIDATION_TYPE
TYPE <<- 'validation_data' 
TRIM_TYPE <<- VALIDATION_TRIM_TYPE
GENE_NAME <<- VALIDATION_GENE_NAME
PRODUCTIVITY <<- VALIDATION_PRODUCTIVITY

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH, '/scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_evaluation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/analysis_scripts/pwm_profile_functions.R'))

validation_data = aggregate_validation_data(directory = VALIDATION_DATA_DIR)
validation_data = get_model_feature_scores(validation_data, pwm)

rel_importance_model = fit_rel_importance_model(validation_data)
loss = evaluate_loss(validation_data, rel_importance_model)
write_rel_importance_result_dt(coef(rel_importance_model), loss, output_file_name)
