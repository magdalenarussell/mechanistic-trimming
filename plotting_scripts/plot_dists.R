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

TRIM_TYPE <<- args[1]
stopifnot(TRIM_TYPE == 'v_trim')

MOTIF_TYPE <<- args[2]
stopifnot(MOTIF_TYPE %in% c('bounded', 'unbounded', 'unbounded_no_pnuc'))

NCPU <<- as.numeric(args[3])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- args[4]
GENE_WEIGHT_TYPE <<- args[5]

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[6])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[7])

UPPER_TRIM_BOUND <<- as.numeric(args[8]) 
LOWER_TRIM_BOUND <<- RIGHT_NUC_MOTIF_COUNT - 2 

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
