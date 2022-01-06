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
# model_types_neat = filter_model_types(remove_types_with_string = c('gc_content'))
model_types_neat = filter_model_types() 
# load evaluation file
file_path = get_model_evaluation_file_name(TYPE)
eval_data = fread(file_path)

# process evaluation file
eval_data = process_model_evaluation_file(eval_data, model_types_neat)
eval_data = eval_data[motif_type == MOTIF_TYPE]

# filter models (for now, comparing only 4x4 motifs and terminal melting 5' length of 10 (or NA))
motifs = eval_data[model_type %like% 'motif' & motif_length_5_end == 4 & motif_length_3_end == 4 & terminal_melting_5_end_length %in% c(NA, 10)]
terminal = eval_data[!(model_type %like% 'motif') & model_type %like% 'terminal' & motif_length_5_end %in% c(0, 4) & terminal_melting_5_end_length %in% c(NA, 10)]
distance = eval_data[model_type == 'distance']

together = rbind(motifs, terminal, distance)

# get total model term count for each model
together = get_term_count(together, 4, 4)

# create a new column with melting_type
# together = get_terminal_melting_calculation_type(together)

plot = ggplot(together) +
    geom_point(aes(y = get(TYPE), x = terms, color = model_type), size = 5)+
    theme_cowplot(font_family = 'Arial') + 
    xlab('Total number of terms') +
    ylab('Conditional log loss') +
    theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

path = get_model_eval_file_path(TYPE)
file_name = paste0(path, '/neat_', TYPE, '_term_count_scatter_no_label.pdf')
ggsave(file_name, plot = plot, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
