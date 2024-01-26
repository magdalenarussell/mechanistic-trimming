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

all = read_frames_data()
all_subset = all[frame_type == 'Out' | frame_stop == TRUE]

cols = colnames(all_subset)[!(colnames(all_subset) %like% 'frame')]
cols = cols[!(cols %like% 'seq_len')]

all_subset = all_subset[, ..cols]

# get oriented full sequences
genes = get_gene_order(GENE_NAME)
whole_nucseq = get_oriented_whole_nucseqs()
seqs = data.table()
for (g in genes){
    gene_seqs = whole_nucseq[substring(gene, 4, 4) %in% toupper(substring(g, 1, 1))]
    colnames(gene_seqs) = c(g, paste0(g, '_sequence'))
    seqs = rbind(seqs, tog, fill = TRUE)
}

v_seqs = seqs[is.na(j_gene), v_ind := 1:.N, by = .(v_gene)][,-c('j_gene', 'j_gene_sequence', 'j_gene')]
v_seqs = v_seqs[v_ind == 1][!is.na(v_gene)][, -c('v_ind', 'j_ind', 'v_gene')]

j_seqs = seqs[is.na(v_gene), j_ind := 1:.N, by = .(j_gene)][,-c('v_gene', 'v_gene_sequence', 'v_gene')]
j_seqs = j_seqs[j_ind == 1][!is.na(j_gene)][, -c('j_ind', 'v_ind', 'j_gene')]

all_subset = merge(all_subset, v_seqs, by = 'v_gene')
all_subset = merge(all_subset, j_seqs, by = 'j_gene')

# get NT context
all_subset$vj_insert = 0
motif_data = get_all_nuc_contexts(all_subset, subject_id=NULL, gene_type = GENE_NAME, trim_type = TRIM_TYPE)

# Filter motif data for possible sites
motif_data = filter_motif_data_for_possible_sites(motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE)

# process data for modeling
processed = inner_aggregation_processing(motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE)

processed = subset_processed_data(processed, trim_type = TRIM_TYPE, gene_type = GENE_NAME)
not_cols = c('weighted_observation', 'count', 'total_tcr')
cols = colnames(processed)[!(colnames(processed) %in% not_cols)]
processed = processed[, ..cols]

filename = processed_data_path()
file_name = str_replace(filename, 'processed_data', 'pre_simulation_processed_data')
fwrite(processed, file_name, sep = '\t')
cat(file_name)
