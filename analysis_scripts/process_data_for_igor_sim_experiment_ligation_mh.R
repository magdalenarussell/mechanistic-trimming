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

ANNOTATION_TYPE <<- 'igor_sim_alpha'
PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- 2
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_ligation-mh'

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# Read raw data
require(vroom)
files = list.files(TCR_REPERTOIRE_DATA_igor_sim_alpha, full.names = TRUE)
all_data = vroom(files)
all_data = data.table(all_data)

all_data = all_data[vj_insert == 0]

# First, generate dataset where V-choice is independent of J-choice
set.seed(55)
vcols = c('seq_index', 'sample_name', 'v_gene', 'v_trim')
jcols = c('seq_index', 'sample_name', 'j_gene', 'j_trim')
v_data = all_data[, ..vcols]
j_data = all_data[, ..jcols]

shuffled_v_data = v_data[sample(nrow(v_data)), ]
shuffled_j_data = j_data[sample(nrow(j_data)), ]
shuffled_v_data[, new_index := seq(.N)]
shuffled_j_data[, new_index := seq(.N)]

together = merge(shuffled_v_data, shuffled_j_data, by = c('new_index', 'sample_name'))
together$productive = FALSE
together$vj_insert = 0 

compiled = compile_data_for_subject(dataset=together, write = FALSE)
compiled_processed = inner_aggregation_processing(compiled, gene_type=GENE_NAME, trim_type=TRIM_TYPE)

ANNOTATION_TYPE <<- 'igor_sim_alpha_independent_vj_choice'
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

compiled_processed = subset_processed_data(compiled_processed, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

filename = processed_data_path()
fwrite(compiled_processed, filename, sep = '\t')


ANNOTATION_TYPE <<- 'igor_sim_alpha'

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# Second, generate dataset where trimming distributions for each V-gene and J-gene are shuffled
vcols = c('v_gene')
jcols = c('j_gene')

v = unique(all_data[, ..vcols])
j = unique(all_data[, ..jcols])

colnames(v) = 'old_v_gene'
colnames(j) = 'old_j_gene'

v[, v_gene := sample(old_v_gene, replace = FALSE)]
j[, j_gene := sample(old_j_gene, replace = FALSE)]

trim_compiled = merge(v, all_data, by.x = 'old_v_gene', by.y = 'v_gene')
trim_compiled = merge(j, trim_compiled, by.x = 'old_j_gene', by.y = 'j_gene')

trim_compiled = compile_data_for_subject(dataset=trim_compiled, write = FALSE)
trim_compiled_processed = inner_aggregation_processing(trim_compiled, gene_type=GENE_NAME, trim_type=TRIM_TYPE)

ANNOTATION_TYPE <<- 'igor_sim_alpha_shuffled_trimming'
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

trim_compiled_processed = subset_processed_data(trim_compiled_processed, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

filename = processed_data_path()
fwrite(trim_compiled_processed, filename, sep = '\t')

ANNOTATION_TYPE <<- 'igor_sim_random_sequence'

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))


rand = compile_data_for_subject(dataset=all_data, write = FALSE)
rand_processed = inner_aggregation_processing(rand, gene_type=GENE_NAME, trim_type=TRIM_TYPE)

rand_processed = subset_processed_data(rand_processed, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

filename = processed_data_path()
fwrite(rand_processed, filename, sep = '\t')

# fourth, generate dataset where trimming distributions for each V-gene and J-gene are shuffled after the trimming adjustment procedure (resulting ligationMH scenarios may not actually be possible given germline)

ANNOTATION_TYPE <<- 'igor_sim_alpha'
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R')) 

filename = processed_data_path()
processed = fread(filename)

ANNOTATION_TYPE <<- 'igor_sim_alpha_shuffled_adjusted_trimming'
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

vcols = c('v_gene')
jcols = c('j_gene')

v = unique(processed[, ..vcols])
j = unique(processed[, ..jcols])

colnames(v) = 'old_v_gene'
colnames(j) = 'old_j_gene'

v[, v_gene := sample(old_v_gene, replace = FALSE)]
j[, j_gene := sample(old_j_gene, replace = FALSE)]

trim_compiled = merge(v, processed, by.x = 'old_v_gene', by.y = 'v_gene')
trim_compiled = merge(j, trim_compiled, by.x = 'old_j_gene', by.y = 'j_gene')

cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh', 'count', 'total_tcr')
trim_compiled_subset = trim_compiled[, ..cols]
trim_compiled_subset$vj_insert = 0

# Get oriented full sequences and group genes by common features
genes = get_gene_order(GENE_NAME)
whole_nucseq = get_oriented_whole_nucseqs()
final = data.table()
for (g in genes){
    gene_seqs = whole_nucseq[substring(gene, 4, 4) %in% toupper(substring(g, 1, 1))]
    colnames(gene_seqs) = c(g, paste0(g, '_sequence'))
    final = rbind(final, tog, fill = TRUE)
}
gcols = c('v_gene', 'j_gene', 'v_gene_sequence', 'j_gene_sequence')
final = unique(final[, ..gcols])

v_groups = final[!is.na(v_gene)]
non_vdup = !duplicated(v_groups$v_gene)
v_groups = v_groups[non_vdup, ]

j_groups = final[!is.na(j_gene)]
non_jdup = !duplicated(j_groups$j_gene)
j_groups = j_groups[non_jdup, ]

trim_compiled_subset = merge(trim_compiled_subset, v_groups[, c('v_gene', 'v_gene_sequence')], by = 'v_gene')
trim_compiled_subset = merge(trim_compiled_subset, j_groups[, c('j_gene', 'j_gene_sequence')], by = 'j_gene')

# Get NT context for each scenario
motif_data = get_all_nuc_contexts(trim_compiled_subset, NULL, gene_type = GENE_NAME, trim_type = TRIM_TYPE)

cols = colnames(trim_compiled_subset)[!(colnames(trim_compiled_subset) %like% 'sequence')]
cols = cols[!(cols %in% c('vj_insert', 'total_tcr'))]

subset = merge(trim_compiled_subset, motif_data, by = cols)

trim_compiled_processed = inner_aggregation_processing(subset, gene_type=GENE_NAME, trim_type=TRIM_TYPE)
trim_compiled_processed = subset_processed_data(trim_compiled_processed, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

filename = processed_data_path()
fwrite(trim_compiled_processed, filename, sep = '\t')

