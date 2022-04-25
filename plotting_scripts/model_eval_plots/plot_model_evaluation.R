source('config/config.R')

library(ggplot2)
library(cowplot)
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

ANNOTATION_TYPE<<- args[1]
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

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[6]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

LOWER_TRIM_BOUND <<- 2
UPPER_TRIM_BOUND <<- args[7] 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- args[8]

TYPE <<- args[9]
stopifnot(TYPE %in% c('log_loss', 'expected_log_loss', 'aic', 'raw_loss', 'old_loss_cv', 'log_loss_j_gene', 'v_gene_family_loss', 'full_v_gene_family_loss'))

source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/model_evaluation_functions.R')

eval_results = compile_evaluation_results(TYPE)

model_types = filter_model_types(remove_types_with_string = c('NN', 'combo', 'distance_terminal_melting', 'motif_terminal_melting', 'gc_content', 'dna_shape'))

model_types = model_types[(model_types %like% 'motif') | (model_types %like% 'shape')]

for (model_type in model_types){
    if (model_type %like% 'two_side' | model_type %like% 'two-side'){
        left_melt = LEFT_SIDE_TERMINAL_MELT_LENGTH
    } else {
        left_melt = NA
    }
    
    plot_model_evaluation_heatmap(eval_results, TYPE, model_type = model_type, terminal_melting_5_end_length_filter= left_melt)
}

