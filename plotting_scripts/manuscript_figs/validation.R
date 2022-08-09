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


ANNOTATION_TYPE <<- 'igor' 
TRIM_TYPE <<- 'v_trim'
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
PRODUCTIVITY <<- 'nonproductive'

MOTIF_TYPE <<- 'unbounded' 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(2)

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- 'p_gene_marginal' 
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

UPPER_TRIM_BOUND <<- as.numeric(14) 
LOWER_TRIM_BOUND <<- 2 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(10)

TYPE <<- 'log_loss'

source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/model_evaluation_functions.R')

all_eval_results= compile_evaluation_results(TYPE)
# get model types
orig_model_types = c("motif_two-side-base-count-beyond", 'null')
# pre-filter data
all_eval_results1 = process_model_evaluation_file(all_eval_results, orig_model_types, left_motif_size_filter = LEFT_NUC_MOTIF_COUNT, right_motif_size_filter = RIGHT_NUC_MOTIF_COUNT, terminal_melting_5_end_length_filter = c(NA, LEFT_SIDE_TERMINAL_MELT_LENGTH))
eval_data_murugan = process_model_evaluation_file(all_eval_results, 'motif', 2, 4, NA)
all_eval_results = rbind(all_eval_results1, eval_data_murugan)
all_eval_results = all_eval_results[motif_type == MOTIF_TYPE]
setnames(all_eval_results, 'log_loss', 'loss')

# annotation_types = c('validation_data_alpha', 'validation_data_beta', 'validation_data_gamma')
annotation_types = c('validation_data_alpha', 'validation_data_beta', 'validation_data_gamma')
trim_types = c('j_trim', 'v_trim')

for (trim_type in trim_types){
    for (type in annotation_types) {
        ANNOTATION_TYPE <<- type
        TRIM_TYPE <<- trim_type
        GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

        source('scripts/model_evaluation_functions.R')
        source('plotting_scripts/plotting_functions.R')
        source('plotting_scripts/model_evaluation_functions.R')

        temp_eval_results = compile_evaluation_results(type)
        setnames(temp_eval_results, type, 'loss')
        temp_eval_results$loss_type = type
        temp_eval_results$trim_type = trim_type
        all_eval_results = rbind(all_eval_results, temp_eval_results, fill = TRUE)
    }
}

# get model types
orig_model_types = c("motif_two-side-base-count-beyond", "motif", 'null')
new_model_types = c('1x2motif + two-side-base-count-beyond (12 params)', '2x4motif (18 params)', 'null (0 params)')

all_eval_results$model_type = mapvalues(all_eval_results$model_type, from = orig_model_types, to = new_model_types) 

neat_names = make_model_names_neat(new_model_types)
colors = set_color_palette(c(neat_names, '2x4motif (18 params)'), with_params = TRUE)

all_eval_results$loss_type = paste0(all_eval_results$loss_type, '_', all_eval_results$trim_type)
loss_types = unique(all_eval_results$loss_type)
nice_loss_types = c('full TCRB V-gene\ntraining dataset', 'TCRA J-gene\ntesting dataset', 'TCRB J-gene\ntesting dataset', 'TCRG J-gene\ntesting dataset', 'TCRA V-gene\ntesting dataset', 'TCRB V-gene\ntesting dataset', 'TCRG V-gene\ntesting dataset')
loss_order = c('full TCRB V-gene\ntraining dataset', 'TCRB V-gene\ntesting dataset', 'TCRA V-gene\ntesting dataset', 'TCRG V-gene\ntesting dataset', 'TCRB J-gene\ntesting dataset', 'TCRA J-gene\ntesting dataset', 'TCRG J-gene\ntesting dataset')
all_eval_results$loss_type = mapvalues(all_eval_results$loss_type, from = loss_types, to=nice_loss_types)

plot = plot_model_evaluation_loss_paracoord(all_eval_results, model_type_list = c(new_model_types, '2x4motif'), pre_filter = TRUE, left_motif_size_filter = LEFT_NUC_MOTIF_COUNT, right_motif_size_filter = RIGHT_NUC_MOTIF_COUNT, terminal_melting_5_end_length_filter = c(NA, LEFT_SIDE_TERMINAL_MELT_LENGTH), loss_bound = c(2.04, 3.01), color_palette = colors, write_plot = FALSE, expand_var = 3.1, loss_order = loss_order) +
    ylab('Expected per-sequence log loss\n')

ANNOTATION_TYPE <<- 'igor'
path = get_manuscript_path()
file_name = paste0(path, '/loss_compare_validation.pdf')
ggsave(file_name, plot = plot, width = 40, height = 20, units = 'in', dpi = 750, device = cairo_pdf)


