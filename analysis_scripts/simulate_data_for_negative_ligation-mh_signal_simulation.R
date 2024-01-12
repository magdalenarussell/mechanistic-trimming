source('mechanistic-trimming/config/config.R')
source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))

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

ANNOTATION_TYPE <<- 'igor_sim_alpha'
PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- 2
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond'

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# get predictions (i.e. probabilities) from no-MH model
probs_path = get_validation_predictions_file_path(L2='False', validation_annotation='pre_simulation')
probs = fread(probs_path)
not_cols = c('count', 'weighted_observation', 'total_tcr', 'v_gene_group_j_gene_group_int', 'v_trim_j_trim_ligation_mh', 'v_trim_j_trim_ligation_mh_int', 'v_gene_group_j_gene_group')
cols = colnames(probs)[!(colnames(probs) %in% not_cols)]
cols = cols[!(cols %like% 'pos1_')]
cols = cols[!(cols %like% 'pos2_')]

probs = probs[, ..cols]

ANNOTATION_TYPE <<- 'no_mh_model_sim_alpha'
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_ligation-mh'
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# total number of sequences
total = 1700000

# get V and J gene sequences
vgenes = unique(probs$v_gene_group)
jgenes = unique(probs$j_gene_group)

# simulate single sequence
set.seed(123) 

sim_data = data.table(v_gene_group = sample(vgenes, total, replace = TRUE),
                      j_gene_group = sample(jgenes, total, replace = TRUE))
sim_data[, row_id := .I]

# Join the data.tables
joined = sim_data[probs, on = .(v_gene_group, j_gene_group), allow.cartesian=TRUE]

sample_trimming_configuration <- function(subset_dt, cols) {
  subset_dt[sample(.N, size = 1, prob = predicted_prob), ..cols]
}

# Apply the function to each group and get the results
sampled = joined[, sample_trimming_configuration(.SD, c('v_trim', 'j_trim', 'ligation_mh')), by = .(v_gene_group, j_gene_group, row_id)]

# get oriented full sequences
genes = get_gene_order(GENE_NAME)
whole_nucseq = get_oriented_whole_nucseqs()
seqs = data.table()
for (g in genes) {
    gene_seqs = whole_nucseq[substring(gene, 4, 4) %in% toupper(substring(g, 1, 1))]
    gene_groups = get_common_genes_from_seqs(gene_seqs, g)
    gene_groups = merge(gene_groups, whole_nucseq, by.x = g, by.y = 'gene')
    seqs = rbind(seqs, gene_groups, fill = TRUE)
}

v_seqs = seqs[is.na(j_gene_group), v_ind := 1:.N, by = .(v_gene_group)][,-c('j_gene', 'j_gene_sequence', 'j_gene_group')]
v_seqs = v_seqs[v_ind == 1][!is.na(v_gene_group)][, -c('v_ind', 'j_ind', 'v_gene')]

j_seqs = seqs[is.na(v_gene_group), j_ind := 1:.N, by = .(j_gene_group)][,-c('v_gene', 'v_gene_sequence', 'v_gene_group')]
j_seqs = j_seqs[j_ind == 1][!is.na(j_gene_group)][, -c('j_ind', 'v_ind', 'j_gene')]

sampled = merge(sampled, v_seqs, by = 'v_gene_group')
sampled = merge(sampled, j_seqs, by = 'j_gene_group')

sampled$vj_insert = 0
motif_data = get_all_nuc_contexts(sampled, subject_id=NULL, gene_type = GENE_NAME, trim_type = TRIM_TYPE)

# Filter motif data for possible sites
motif_data = filter_motif_data_for_possible_sites(motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE)

# process data for modeling
processed = inner_aggregation_processing(motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE)

processed = subset_processed_data(processed, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

filename = processed_data_path()
fwrite(processed, filename, sep = '\t')
