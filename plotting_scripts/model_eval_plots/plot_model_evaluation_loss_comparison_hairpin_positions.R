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

MOTIF_TYPE <<- 'unbounded' 

NCPU <<- as.numeric(args[4])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[5]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[6])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[7])

UPPER_TRIM_BOUND <<- as.numeric(args[8]) 
LOWER_TRIM_BOUND <<- 2 

MODEL_TYPE <<- args[9]
stopifnot(MODEL_TYPE != 'null')

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[10])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

TYPE <<- 'log_loss'
all_types = c('log_loss', 'expected_log_loss', 'v_gene_family_loss', 'log_loss_j_gene', 'full_v_gene_family_loss')

source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/model_evaluation_functions.R')

all_eval_results = data.table()
for (type in all_types) {
    temp_eval_results = compile_evaluation_results(type)
    setnames(temp_eval_results, type, 'loss')
    temp_eval_results$loss_type = type
    all_eval_results = rbind(all_eval_results, temp_eval_results, fill = TRUE)
}

all_eval_results[loss_type %like% 'v_gene_family_loss', loss_type := paste0(loss_type, ', cluster ', held_out_clusters)]

motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
nice_hairpin_names = make_hairpin_names_neat(motif_types) 
colors = set_color_palette(c(nice_hairpin_names))

plot_model_evaluation_loss_paracoord(all_eval_results, model_type_list = MODEL_TYPE, left_motif_size_filter = LEFT_NUC_MOTIF_COUNT, right_motif_size_filter = RIGHT_NUC_MOTIF_COUNT, terminal_melting_5_end_length_filter = c(NA, LEFT_SIDE_TERMINAL_MELT_LENGTH), loss_bound = c(1.82, 2.35), color_palette = colors, same_motif_type = FALSE, custom_name = paste0(MODEL_TYPE, '_hairpin_nick_analysis'), plot_size = c(40, 25))
