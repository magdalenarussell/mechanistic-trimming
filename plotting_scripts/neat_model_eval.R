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
stopifnot(TRIM_TYPE == 'v_trim')

MOTIF_TYPE <<- args[3] 
stopifnot(MOTIF_TYPE %in% c('bounded', 'unbounded', 'unbounded_no_pnuc'))

NCPU <<- as.numeric(args[4])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[5]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

LOWER_TRIM_BOUND <<- 2
UPPER_TRIM_BOUND <<- args[6] 

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

motifs = eval_data[model_type %like% 'motif' & motif_length_5_end == 4 & motif_length_3_end == 4]
terminal = eval_data[!(model_type %like% 'motif') & model_type %like% 'terminal' & motif_length_5_end == 4]
distance = eval_data[model_type == 'distance']

together = rbind(motifs, terminal, distance)

together[, terms := 0]
together[model_type %like% 'motif', terms := terms + 8*4]
together[model_type %like% 'distance', terms := terms + 1]
together[model_type %like% 'terminal_melting', terms := terms + 1]
together$model_type = factor(together$model_type, levels = c('motif', 'terminal_melting', 'distance', 'distance_terminal_melting', 'motif_terminal_melting', 'motif_distance', 'motif_distance_terminal_melting'))

together$neat_model_type = mapvalues(together$model_type, from = c('motif', 'terminal_melting', 'distance', 'distance_terminal_melting', 'motif_terminal_melting', 'motif_distance', 'motif_distance_terminal_melting'), to = c('motif', 'sequence breathing', 'distance', 'distance + sequence breathing', 'motif + sequence breathing', 'motif + distance', 'motif + distance + sequence breathing'))

plot = ggplot(together) +
    geom_text_repel(aes(y = log_loss, x = terms, color = model_type, label = neat_model_type), size = 7, direction = 'y', hjust = 0) + 
    theme_cowplot(font_family = 'Arial') + 
    xlab('Total number of terms') +
    ylab('Conditional log loss') +
    theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank(), legend.position = 'none') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_x_continuous(breaks = seq(1, 60, 2), limits = c(0, 60))

path = get_model_eval_file_path('log_loss')
file_name = paste0(path, '/neat_log_loss_term_count_scatter.pdf')
ggsave(file_name, plot = plot, width = 14, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
