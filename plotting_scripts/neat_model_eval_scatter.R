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

TYPE <<- args[7]

source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')

model_type_files = list.files(path = 'scripts/model_formula_functions/')
model_types = str_sub(model_type_files[model_type_files != '_ignore'], end = -3)
model_types_neat = model_types[!grepl('gc_content', model_types)]
model_types_neat = model_types_neat[order(model_types_neat)]

file_path = get_model_evaluation_file_name(TYPE)
eval_data = fread(file_path)
eval_data = eval_data[model_type %in% model_types_neat]

motifs = eval_data[model_type %like% 'motif' & motif_length_5_end == 4 & motif_length_3_end == 4 & terminal_melting_5_end_length %in% c(NA, 10)]
terminal = eval_data[!(model_type %like% 'motif') & model_type %like% 'terminal' & motif_length_5_end %in% c(0, 4) & terminal_melting_5_end_length %in% c(NA, 10)]
distance = eval_data[model_type == 'distance']

together = rbind(motifs, terminal, distance)

together[, terms := 0]
together[model_type %like% 'motif', terms := terms + 8*4]
together[model_type %like% 'distance', terms := terms + 16]
together[model_type %like% 'terminal_melting', terms := terms + 1]
together[model_type %like% 'two_side_terminal_melting', terms := terms + 1]

together[model_type %like% 'terminal_melting_', melting_type := sapply(model_type, function(x) tail(str_split(x, '_')[[1]], 1))]
together[model_type %like% 'terminal_melting' & is.na(melting_type), melting_type := 'simple']
together[is.na(melting_type), melting_type := 'NA']

together[model_type %like% 'terminal_melting_NN', model_type := substring(model_type, 1, nchar(model_type)- 3)]
together[model_type %like% 'terminal_melting_combo', model_type := substring(model_type, 1, nchar(model_type)- 6)]

plot = ggplot(together) +
    geom_point(aes(y = get(TYPE), x = terms, color = model_type, shape = melting_type), size = 5)+
    theme_cowplot(font_family = 'Arial') + 
    xlab('Total number of terms') +
    ylab('Conditional log loss') +
    theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

path = get_model_eval_file_path(TYPE)
file_name = paste0(path, '/neat_', TYPE, '_term_count_scatter_no_label.pdf')
ggsave(file_name, plot = plot, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
