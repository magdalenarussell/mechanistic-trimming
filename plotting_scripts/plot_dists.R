source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(ggplot2)
library(cowplot)
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

MODEL_GROUP <<- args[6]
GENE_WEIGHT_TYPE <<- args[7]

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[8])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[9])

UPPER_TRIM_BOUND <<- as.numeric(args[10]) 

MODEL_TYPE <<- args[11]

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[12])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')

# Read in dist data
predicted_trims = get_predicted_distribution_data() 

for (gene in unique(predicted_trims$gene)){
    plot_predicted_trimming_dists(predicted_trims, gene)
    plot_model_residual_boxplot(predicted_trims, gene)
}

plot_all_model_residuals(predicted_trims)
