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
trim_types = list.files(path = 'scripts/gene_specific_functions/')
trim_types = str_sub(trim_types, end = -3)
stopifnot(TRIM_TYPE %in% trim_types)

PRODUCTIVITY <<- 'nonproductive'

MOTIF_TYPE <<- 'unbounded' 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(2)

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

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

# get model types
orig_model_types = c("motif_distance", "motif_two-side-base-count-beyond", "null", "motif", "distance", "two-side-base-count", 'dna_shape-std', 'linear-distance')
new_model_types = c('1x2motif + categorical-distance', '1x2motif + two-side-base-count-beyond', 'null', "1x2motif", "categorical-distance", 'two-side-base-count', '1x2DNA-shape', 'linear-distance')

# pre-filter data
subset_eval_data = process_model_evaluation_file(all_eval_results, orig_model_types, left_motif_size_filter = LEFT_NUC_MOTIF_COUNT, right_motif_size_filter = RIGHT_NUC_MOTIF_COUNT, terminal_melting_5_end_length_filter = c(NA, LEFT_SIDE_TERMINAL_MELT_LENGTH))
    
subset_eval_data$model_type = mapvalues(subset_eval_data$model_type, from = orig_model_types, to = new_model_types) 

neat_names = make_model_names_neat(new_model_types)
colors = set_color_palette(c(neat_names, '2x4motif'))

eval_data_murugan = process_model_evaluation_file(all_eval_results, 'motif', 2, 4, NA)
eval_data_murugan$model_type = mapvalues(eval_data_murugan$model_type, from = 'motif', to = '2x4motif')

eval_tog = rbind(subset_eval_data, eval_data_murugan)

#filter out extra "most different" held-out group trials
eval_tog = eval_tog[!(loss_type == "v_gene_family_loss, cluster 2" | loss_type == "full_v_gene_family_loss, cluster 2")]
eval_tog = eval_tog[!(loss_type == "v_gene_family_loss, cluster 3" | loss_type == "full_v_gene_family_loss, cluster 3")]
eval_tog = eval_tog[!(loss_type == "v_gene_family_loss, cluster 4" | loss_type == "full_v_gene_family_loss, cluster 4")]


loss_types = unique(eval_tog$loss_type)
nice_loss_types = c('full V-gene\ntraining\ndataset', 'many held-out\nsubsets of\nV-gene\ntraining\ndataset', '\"most different\"\ncluster of\nV-genes\n(terminal seqs)', 'full J-gene\ndataset', '\"most different\"\ncluster of\nV-genes\n(full seqs)')

eval_tog$loss_type = mapvalues(eval_tog$loss_type, from = loss_types, to=nice_loss_types)

plot = plot_model_evaluation_loss_paracoord(eval_tog, model_type_list = c(new_model_types, '2x4motif'), pre_filter = TRUE, left_motif_size_filter = LEFT_NUC_MOTIF_COUNT, right_motif_size_filter = RIGHT_NUC_MOTIF_COUNT, terminal_melting_5_end_length_filter = c(NA, LEFT_SIDE_TERMINAL_MELT_LENGTH), loss_bound = c(1.97, 2.7), color_palette = colors, write_plot = FALSE, expand_var = 2.3) +
    ylab('Expected per-sequence log loss\n')

path = get_manuscript_path()
file_name = paste0(path, '/loss_compare.pdf')
ggsave(file_name, plot = plot, width = 32, height = 20, units = 'in', dpi = 750, device = cairo_pdf)


