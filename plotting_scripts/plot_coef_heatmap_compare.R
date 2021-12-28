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

MOTIF_TYPE <<- args[3] 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[4])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- 'all_subjects' 
GENE_WEIGHT_TYPE <<- args[5]

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[6])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[7])

UPPER_TRIM_BOUND <<- as.numeric(args[8]) 
LOWER_TRIM_BOUND <<- RIGHT_NUC_MOTIF_COUNT - 2 

MODEL_TYPE <<- args[9]
model_type1 = MODEL_TYPE

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- args[10]

if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(LEFT_SIDE_TERMINAL_MELT_LENGTH)
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

left_side_terminal_melt_length1 = LEFT_SIDE_TERMINAL_MELT_LENGTH

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')

# Read in model coefficient data 
pwm1 = get_model_coefficient_data() 

# get model type for the second heatmap 
MODEL_TYPE <<- args[11]
model_type2 = MODEL_TYPE

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- args[12]

if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(LEFT_SIDE_TERMINAL_MELT_LENGTH)
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

left_side_terminal_melt_length2 = LEFT_SIDE_TERMINAL_MELT_LENGTH

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')

# Read in model coefficient data 
pwm2= get_model_coefficient_data() 

file_name = get_model_coef_heatmap_compare_file_name(model_type1, model_type2, left_side_terminal_melt_length1, left_side_terminal_melt_length2)
plot_model_coefficient_heatmap_compare(pwm1, pwm2, file_name, with_values = TRUE ,limits = NULL, write_plot = TRUE)
