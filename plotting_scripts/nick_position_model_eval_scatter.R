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

NCPU <<- as.numeric(args[4])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- 'all_subjects'
MOTIF_TYPE <<- 'unbounded'

GENE_WEIGHT_TYPE <<- args[5]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

LOWER_TRIM_BOUND <<- 2
UPPER_TRIM_BOUND <<- as.numeric(args[6]) 

TYPE <<- args[7]

source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/model_evaluation_functions.R')

# get model types
model_types_neat = filter_model_types(remove_types_with_string = c('gc_content', 'NN', 'combo'))

# load evaluation file
file_path = get_model_evaluation_file_name(TYPE, intermediate = FALSE)
eval_data = fread(file_path)

# process evaluation file
together = process_model_evaluation_file(eval_data, model_types_neat, left_motif_size_filter = 4, right_motif_size_filter = 4, terminal_melting_5_end_length_filter = c(NA, 10))

# get total model term count for each model
together = get_term_count(together, 4, 4)

# get nick positions 
together = get_pnuc_count(together)
together = together[order(pnuc_count)]
together$nick_position = factor(together$nick_position, levels = list(unique(together$nick_position), unique(together$pnuc_count))[[1]])

plot = ggplot(together) +
    geom_line(aes(y = get(TYPE), x = nick_position, color = model_type, group = model_type), size = 2) +
    geom_point(aes(y = get(TYPE), x = nick_position, color = model_type), size = 5) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Nick position') +
    ylab('Conditional log loss') +
    theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

path = get_model_eval_file_path(TYPE)
file_name = paste0(path, '/neat_', TYPE, '_nick_position_scatter.pdf')
ggsave(file_name, plot = plot, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
