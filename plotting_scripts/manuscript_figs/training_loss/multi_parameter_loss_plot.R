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

LOSS_GENE_WEIGHT <<- 'p_gene_given_subject'

source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/model_evaluation_functions.R')

# get all loss results
all_eval_results = data.table()
for (type in all_types) {
    temp_eval_results = compile_evaluation_results(type)
    setnames(temp_eval_results, type, 'loss')
    temp_eval_results$loss_type = type
    all_eval_results = rbind(all_eval_results, temp_eval_results, fill = TRUE)
}

# assign cluster for "most different" sequence protocol
all_eval_results[loss_type %like% 'v_gene_family_loss', loss_type := paste0(loss_type, ', cluster ', held_out_clusters)]

# get model types, and make them neat
orig_model_types = c("motif_two-side-base-count-beyond", "null", "motif", "two-side-base-count", 'dna_shape-std', 'linear-distance', 'motif_linear-distance')
new_model_types = c('1x2motif + two-side\nbase-count beyond\n(12 params)', 'null (0 params)', "1x2motif (9 params)", 'two-side base-count\n(3 params)', '1x2DNA-shape\n(13 params)', 'distance (1 param)', '1x2motif + distance\n(10 params)')

# pre-filter data
subset_eval_data = process_model_evaluation_file(all_eval_results, orig_model_types, left_motif_size_filter = LEFT_NUC_MOTIF_COUNT, right_motif_size_filter = RIGHT_NUC_MOTIF_COUNT, terminal_melting_5_end_length_filter = c(NA, LEFT_SIDE_TERMINAL_MELT_LENGTH))
    
# reassign model names
subset_eval_data$model_type = mapvalues(subset_eval_data$model_type, from = orig_model_types, to = new_model_types) 

# make model names neater and set palette colors
neat_names = make_model_names_neat(new_model_types)
colors = set_color_palette(c(neat_names, '2x4motif (18 params)'), with_params = TRUE)

# get model evaluation results for the 2x4 motif model
eval_data_murugan = process_model_evaluation_file(all_eval_results, 'motif', 2, 4, NA)
eval_data_murugan$model_type = mapvalues(eval_data_murugan$model_type, from = 'motif', to = '2x4motif (18 params)')

# combine data sets
eval_tog = rbind(subset_eval_data, eval_data_murugan)

# filter out model evaluation types
eval_tog = eval_tog[!(loss_type == "v_gene_family_loss, cluster 2" | loss_type == "full_v_gene_family_loss, cluster 2")]
eval_tog = eval_tog[!(loss_type == "v_gene_family_loss, cluster 3" | loss_type == "full_v_gene_family_loss, cluster 3")]
eval_tog = eval_tog[!(loss_type == "v_gene_family_loss, cluster 4" | loss_type == "full_v_gene_family_loss, cluster 4")]

# order and reassign loss types (neater version)
loss_types = unique(eval_tog$loss_type)
nice_loss_types = c('full V-gene\ntraining\ndataset', 'many held-out\nsubsets of\nV-gene\ntraining\ndataset', '\"most different\"\ncluster of\nV-genes\n(terminal seqs)', 'full J-gene\ndataset', '\"most different\"\ncluster of\nV-genes\n(full seqs)')
eval_tog$loss_type = mapvalues(eval_tog$loss_type, from = loss_types, to=nice_loss_types)

eval_tog = eval_tog[motif_type == MOTIF_TYPE]
cols = c('model_type', 'loss_type', 'loss')
loss_data = eval_tog[, ..cols]
colnames(loss_data) = c('model_type', 'loss_type', 'log_loss')
fwrite(loss_data, paste0(PROJECT_PATH, '/plotting_scripts/manuscript_figs/training_loss/loss.tsv'), sep = '\t')


# create plot
plot = plot_model_evaluation_loss_paracoord(eval_tog, model_type_list = c(new_model_types, '2x4motif'), pre_filter = TRUE, left_motif_size_filter = LEFT_NUC_MOTIF_COUNT, right_motif_size_filter = RIGHT_NUC_MOTIF_COUNT, terminal_melting_5_end_length_filter = c(NA, LEFT_SIDE_TERMINAL_MELT_LENGTH), loss_bound = c(1.98, 2.69), color_palette = colors, write_plot = FALSE, expand_var = 1.2) +
    ylab('Expected per-sequence log loss\n')

# save plot
path = get_manuscript_path()
file_name = paste0(path, '/loss_compare.pdf')
ggsave(file_name, plot = plot, width = 32, height = 17, units = 'in', dpi = 750, device = cairo_pdf)
