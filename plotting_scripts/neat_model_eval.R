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

PRODUCTIVITY <<- args[4]

MOTIF_TYPE <<- args[5] 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[6])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[7]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

LOWER_TRIM_BOUND <<- 2
UPPER_TRIM_BOUND <<- as.numeric(args[8])

source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')

model_type_files = list.files(path = 'scripts/model_formula_functions/')
model_types = str_sub(model_type_files[model_type_files != '_ignore'], end = -3)
model_types_neat = model_types[!grepl('gc_content', model_types)]
model_types_neat = model_types_neat[order(model_types_neat)]

type = 'log_loss'

file_path = get_model_evaluation_file_name(type)
eval_data = fread(file_path)
eval_data = eval_data[model_type %in% model_types_neat]
eval_data = unique(eval_data)
eval_data = eval_data[motif_type == MOTIF_TYPE]
eval_data = eval_data[gene_weight_type == GENE_WEIGHT_TYPE]
eval_data = eval_data[lower_bound == LOWER_TRIM_BOUND]
eval_data = eval_data[upper_bound == UPPER_TRIM_BOUND]


motifs = eval_data[model_type %like% 'motif' & motif_length_5_end == 4 & motif_length_3_end == 4 & terminal_melting_5_end_length %in% c(NA, 10)]
terminal = eval_data[!(model_type %like% 'motif') & model_type %like% 'terminal' & motif_length_5_end %in% c(0, 4) & terminal_melting_5_end_length %in% c(NA, 10)]
distance = eval_data[model_type == 'distance']

together = rbind(motifs, terminal, distance)

together[, terms := 0]
together[model_type %like% 'motif', terms := terms + 8*4]
together[model_type %like% 'distance', terms := terms + 16]
together[model_type %like% 'terminal_melting', terms := terms + 1]
together[model_type %like% 'two_side_terminal_melting', terms := terms + 1]

together$neat_model_type = mapvalues(together$model_type, from = c('motif', 'terminal_melting', 'distance', 'distance_terminal_melting', 
                                                                   'motif_terminal_melting', 'motif_distance', 'motif_distance_terminal_melting', 
                                                                   'two_side_terminal_melting', 'distance_two_side_terminal_melting', 'motif_distance_two_side_terminal_melting', 'motif_two_side_terminal_melting',
                                                                   'terminal_melting_NN', 'distance_terminal_melting_NN', 'motif_terminal_melting_NN', 'motif_distance_terminal_melting_NN', 
                                                                   'two_side_terminal_melting_NN', 'distance_two_side_terminal_melting_NN', 'motif_distance_two_side_terminal_melting_NN', 'motif_two_side_terminal_melting_NN',
                                                                   'terminal_melting_combo', 'distance_terminal_melting_combo', 'motif_terminal_melting_combo', 'motif_distance_terminal_melting_combo', 
                                                                   'two_side_terminal_melting_combo', 'distance_two_side_terminal_melting_combo', 'motif_distance_two_side_terminal_melting_combo', 'motif_two_side_terminal_melting_combo'),
                                     to = c('motif', 'sequence breathing', 'distance', 'distance + sequence breathing', 
                                            'motif + sequence breathing', 'motif + distance', 'motif + distance + sequence breathing', 
                                            'two-side sequence breathing', 'distance + two-side sequence breathing', 'motif + distance + two-side sequence breathing', 'motif + two-side sequence breathing',
                                            'sequence breathing (NN)', 'distance + sequence breathing (NN)', 'motif + sequence breathing (NN)', 'motif + distance + sequence breathing (NN)', 
                                            'two-side sequence breathing (NN)', 'distance + two-side sequence breathing (NN)', 'motif + distance + two-side sequence breathing (NN)', 'motif + two-side sequence breathing (NN)',
                                            'sequence breathing (combo)', 'distance + sequence breathing (combo)', 'motif + sequence breathing (combo)', 'motif + distance + sequence breathing (combo)', 
                                            'two-side sequence breathing (combo)', 'distance + two-side sequence breathing (combo)', 'motif + distance + two-side sequence breathing (combo)', 'motif + two-side sequence breathing (combo)'))


plot = ggplot(together) +
    # geom_point(aes(y = log_loss, x = terms, color = model_type), size = 5)+
    geom_text_repel(aes(y = log_loss, x = terms, color = model_type, label = neat_model_type), size = 5, direction = 'y', hjust = 0) + 
    theme_cowplot(font_family = 'Arial') + 
    xlab('Total number of terms') +
    ylab('Conditional log loss') +
    theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank(), legend.position = 'none') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_x_continuous(breaks = seq(0, 80, 2), limits = c(0, 80))+
    ylim(410000, 500000)

path = get_model_eval_file_path('log_loss')
file_name = paste0(path, '/neat_log_loss_term_count_scatter.pdf')
ggsave(file_name, plot = plot, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
