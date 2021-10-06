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

TRIM_TYPE <<- args[1]
stopifnot(TRIM_TYPE == 'v_trim')

MOTIF_TYPE <<- args[2] 
stopifnot(MOTIF_TYPE %in% c('bounded', 'unbounded', 'unbounded_no_pnuc'))

NCPU <<- as.numeric(args[3])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[4]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

LOWER_TRIM_BOUND <<- 2
UPPER_TRIM_BOUND <<- args[5] 

source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')

model_type_files = list.files(path = 'scripts/model_formula_functions/')
model_types = str_sub(model_type_files[model_type_files != '_ignore'], end = -3)

for (type in c('per_gene', 'log_loss', 'per_gene_per_trim')){
    for (model_type in model_types){
        plot_model_evaluation_heatmap(type, model_type_filter = model_type)
        plot_model_evaluation_scatter(type, model_type_filter = model_type)
    }
    plot_model_evaluation_heatmap(type)
    plot_model_evaluation_scatter(type)
}
