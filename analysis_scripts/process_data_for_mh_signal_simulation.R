source('mechanistic-trimming/config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

ANNOTATION_TYPE <<- 'igor_alpha'
PARAM_GROUP <<- 'nonproductive_v-j_trim'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- 2
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_mh'

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))

# Read processed data
filename = processed_data_path()
motif_data = fread(filename)

# simulate data to replicate what an expected MH signal would yield
PARAM_GROUP <<- 'nonproductive_v-j_trim_mh_sim'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

# get MH prop quantiles
motif_data[, mh_prop_mid_sum := mh_prop_mid_overlap_1 + mh_prop_mid_overlap_2 + mh_prop_mid_overlap_3 + mh_prop_mid_overlap_4]
motif_data[, mh_prop_mid_quantile := cut(mh_prop_mid_sum, quantile(mh_prop_mid_sum, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), labels = c("Q1", "Q2", "Q3", "Q4"))]

motif_data[, mh_prop_outer_sum := mh_prop_up_overlap_1 + mh_prop_up_overlap_2 + mh_prop_up_overlap_3 + mh_prop_up_overlap_4 + mh_prop_down_overlap_1 + mh_prop_down_overlap_2 + mh_prop_down_overlap_3 + mh_prop_down_overlap_4]
motif_data[, mh_prop_outer_quantile := cut(mh_prop_outer_sum, quantile(mh_prop_outer_sum, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE), labels = c("Q1", "Q2", "Q3", "Q4"))]

# define conditions
motif_data[, very_high_mh_condition := (mh_prop_mid_quantile %in% c('Q4')) & (motif_data$mh_prop_outer_quantile %in% c('Q1'))]

motif_data[, high_mh_condition := (mh_prop_mid_quantile %in% c('Q3', 'Q4')) & (motif_data$mh_prop_outer_quantile %in% c('Q1', 'Q2') & very_high_mh_condition == FALSE)]

motif_data[, low_mh_condition := (mh_prop_mid_quantile %in% c('Q1', 'Q2')) & (motif_data$mh_prop_outer_quantile %in% c('Q3', 'Q4') & very_high_mh_condition == FALSE & high_mh_condition == FALSE)]

motif_data[, other_mh_condition := !(very_high_mh_condition == TRUE | high_mh_condition == TRUE | low_mh_condition == TRUE)]

# assign sample probabilities based on these conditions
motif_data[, mh_prob := 0]
motif_data[very_high_mh_condition == TRUE, mh_prob := 1]
motif_data[high_mh_condition == TRUE, mh_prob := 0.8]
motif_data[other_mh_condition == TRUE, mh_prob := 0.1]
motif_data[low_mh_condition == TRUE, mh_prob := 0.05]


# resample counts based on conditions
genes = get_gene_order(GENE_NAME)
trims = get_trim_order(TRIM_TYPE)

# sample size is the number of gene combos
motif_data$cluster = interaction(motif_data[[paste0(genes[1], '_group')]], motif_data[[paste0(genes[2], '_group')]])

# sample proportion of sequences for each individual
size = unique(motif_data$total_tcr)
motif_data[, subsample_total_tcr := size]
vars = c(colnames(motif_data)[colnames(motif_data) %like% 'mh_prop'],
         colnames(motif_data)[colnames(motif_data) %like% 'base_count'],
         colnames(motif_data)[colnames(motif_data) %like% 'motif'], 
         colnames(motif_data)[colnames(motif_data) %like% 'condition'])

cols = c(paste0(genes, '_group'), trims, vars, 'count', 'subsample_total_tcr', 'cluster')
motif_data[, row := seq(1, .N)]
subset = motif_data[motif_data[, sample(.I, size, replace = TRUE, prob = mh_prob)]]
subset[, count := .N, by = .(row)]
subset_final = unique(subset[, ..cols])

# fill in unobserved seq cases
subset_orig_small = motif_data[, ..cols][, -c('count')]
subset_final_small = subset_final[, ..cols][, -c('count')]
unsampled = fsetdiff(subset_orig_small, subset_final_small) 
unsampled$count = 0
subset_final = rbind(subset_final, unsampled)

sample_data = calculate_subject_gene_weight(subset_final, gene_type = GENE_NAME, trim_type = TRIM_TYPE)

filename = processed_data_path()
fwrite(sample_data, filename, sep = '\t')
