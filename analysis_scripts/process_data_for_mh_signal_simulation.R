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
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_interior-mh-count'

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# Read processed data
filename = processed_data_path()
motif_data = fread(filename)

# simulate data to replicate what an expected MH signal would yield
PARAM_GROUP <<- 'nonproductive_v-j_trim_mh_sim'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

# get MH count quantiles
motif_data[, mh_count_mid_sum := mh_count_mid_overlap_1 + mh_count_mid_overlap_2 + mh_count_mid_overlap_3 + mh_count_mid_overlap_4]

# assign sample probabilities based on these conditions
a = 2
tot = sum(log(1+ seq(max(motif_data$mh_count_mid_sum), min(motif_data$mh_count_mid_sum))) + a)
motif_data[, mh_prob := (log(1 + mh_count_mid_sum) + a)/tot]

# resample counts based on conditions
genes = get_gene_order(GENE_NAME)
trims = get_trim_order(TRIM_TYPE)

# sample size is the number of gene combos
motif_data$cluster = interaction(motif_data[[paste0(genes[1], '_group')]], motif_data[[paste0(genes[2], '_group')]])

# sample proportion of sequences for each individual
size = unique(motif_data$total_tcr)
motif_data[, subsample_total_tcr := size]
vars = c(colnames(motif_data)[colnames(motif_data) %like% 'mh_count'],
         colnames(motif_data)[colnames(motif_data) %like% 'base_count'],
         colnames(motif_data)[colnames(motif_data) %like% 'motif'], 
         colnames(motif_data)[colnames(motif_data) %like% 'condition'])

cols = c(paste0(genes, '_group'), trims, vars, 'count', 'subsample_total_tcr', 'cluster')
motif_data[, gene_pair_count := sum(count), by = cluster]
motif_data[, row := seq(1, .N)]
subset = motif_data[motif_data[, .(sampled_rows = sample(.I, size = gene_pair_count[1L], replace = TRUE, prob = mh_prob)), by = cluster]$sampled_rows]
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
