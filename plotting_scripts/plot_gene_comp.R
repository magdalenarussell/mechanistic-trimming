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
stopifnot(TRIM_TYPE == 'v_trim')

MOTIF_TYPE <<- args[3] 
stopifnot(MOTIF_TYPE %in% c('bounded', 'unbounded', 'unbounded_no_pnuc'))

NCPU <<- as.numeric(args[4])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- args[5]
GENE_WEIGHT_TYPE <<- args[6]

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[7])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[8])

UPPER_TRIM_BOUND <<- as.numeric(args[9]) 
LOWER_TRIM_BOUND <<- RIGHT_NUC_MOTIF_COUNT - 2 

MODEL_TYPE <<- args[10]

if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[11])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}


source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')

# Compile data for all subjects
motif_data = aggregate_all_subject_data()

plot_gene_composition(motif_data, weighting = 'uniform')
plot_gene_composition(motif_data, weighting = 'p_gene_marginal')
plot_gene_composition(motif_data, weighting = 'weighted_observation')
plot_gene_composition(motif_data, weighting = 'raw_count')
