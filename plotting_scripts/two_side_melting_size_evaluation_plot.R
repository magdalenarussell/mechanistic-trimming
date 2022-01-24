source('config/config.R')

library(cli)
library(devtools)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(foreach)
library(doParallel)
library(tidyverse)
library(plyr)
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
UPPER_TRIM_BOUND <<- as.numeric(args[7])

TYPE <<- args[8]

source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/model_evaluation_functions.R')

# get model types
model_types_neat = filter_model_types(remove_types_with_string = c('gc_content'))

# load evaluation file
file_path = get_model_evaluation_file_name(TYPE, intermediate = FALSE)
eval_data = fread(file_path)

# process evaluation file
eval_data = process_model_evaluation_file(eval_data, model_types_neat, terminal_melting_5_end_length_filter = c(5, 10, 15, 20), left_motif_size_filter = 4, right_motif_size_filter = 4)
eval_data = eval_data[motif_type == MOTIF_TYPE]
together = eval_data[model_type %like% 'two_side_terminal_melting']

# create a new column with melting_type
# together = get_terminal_melting_calculation_type(together)

plot = ggplot(together) +
    geom_line(aes(y = get(TYPE), x = terminal_melting_5_end_length, color = model_type, group = model_type), size = 2) +
    geom_point(aes(y = get(TYPE), x = terminal_melting_5_end_length, color = model_type), size = 5)+
    theme_cowplot(font_family = 'Arial') + 
    xlab('Number of nucleotides used to calculate 5\' melting') +
    ylab('Expected conditional log loss') +
    theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

path = get_model_eval_file_path(TYPE)
file_name = paste0(path, '/neat_', TYPE, '_two_side_terminal_melting_size_evaluation.pdf')
ggsave(file_name, plot = plot, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
