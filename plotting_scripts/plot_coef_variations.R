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

MOTIF_TYPE <<- 'bounded'
stopifnot(MOTIF_TYPE %in% c('bounded', 'unbounded'))

NCPU <<- as.numeric(args[2])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- args[3]
GENE_WEIGHT_TYPE <<- args[4]

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[5])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[6])

UPPER_TRIM_BOUND <<- 18
LOWER_TRIM_BOUND <<- RIGHT_NUC_MOTIF_COUNT - 2 

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')

# Read in model coefficient data for all subjects in model group, individually 
pwm = get_model_coefficient_data() 

plot_model_coefficient_variations(pwm)
