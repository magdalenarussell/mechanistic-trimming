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

# NOTE: This method is only applicable for models fit across all subjects!
MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[6]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[7])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[8])

UPPER_TRIM_BOUND <<- as.numeric(args[9]) 

MODEL_TYPE <<- 'motif_distance_two_side_terminal_melting' 

if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[10])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('analysis_scripts/bootstrap_analysis_functions.R')

# Compile data for all subjects
motif_data = aggregate_all_subject_data()

# get coefficient count
positions = get_positions()
coef_count = length(positions)*3 + (UPPER_TRIM_BOUND - LOWER_TRIM_BOUND) + 2 

bootstrap = 1
while (bootstrap < 100) {
    boot_data = generate_bootstrap_sample(motif_data, sample_size = length(unique(motif_data$gene))) 
    print(paste0('generated bootstrap sample ', bootstrap, ' of 100'))
    unique = length(unique(boot_data$gene))
    # if (unique < 41){
    #     print(paste0('only ', unique, ' unique genes. Re-generating bootstrap sample.'))
    #     next
    # }
    print(paste0('starting bootstrap fit ', bootstrap, ' of 100'))
    fit_model_bootstrap_genes(boot_data, bootstrap)
    print(paste0('FINISHED bootstrap fit ', bootstrap, ' of 100'))
    bootstrap = bootstrap + 1
}
