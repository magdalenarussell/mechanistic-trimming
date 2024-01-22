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
library(vroom)

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

args = commandArgs(trailingOnly=TRUE)

# get igor parameters
raw_igor_param_path = args[1]
igor_params = fread(raw_igor_param_path)

# get igor simulation data
raw_igor_sim_files = list.files(TCR_REPERTOIRE_DATA_igor_sim_alpha, full.names = TRUE)
igor_sim_data = as.data.table(vroom(raw_igor_sim_files))
igor_sim_data_subset = igor_sim_data[vj_insert == 0]

#########################################
## ADJUSTED IGOR DATA PROCESSING ########
#########################################

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 0 
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 0 
MODEL_TYPE <<- 'adjusted_igor_ligation-mh'
filename = processed_data_path()

# Reformat and filter data
temp_data = reformat_data(igor_sim_data_subset)

# Filter data by productivity and other adaptive data factors
temp_data = filter_by_productivity(temp_data)    
temp_data = adaptive_data_filtering(temp_data)
temp_data = temp_data[v_trim <= UPPER_TRIM_BOUND & v_trim >= LOWER_TRIM_BOUND]
temp_data = temp_data[j_trim <= UPPER_TRIM_BOUND & j_trim >= LOWER_TRIM_BOUND]

# condense data slightly
temp_data_condensed = temp_data[, .N, by = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'productive', 'vj_insert', 'sample_name')]
setnames(temp_data_condensed, 'N', 'count')

# fill in missing scenarios
full = unique(temp_data_condensed[, c('v_gene', 'j_gene')])
full = full[, c(.SD, list(v_trim = LOWER_TRIM_BOUND:UPPER_TRIM_BOUND)), by = .(v_gene, j_gene)]
full = full[, c(.SD, list(j_trim = LOWER_TRIM_BOUND:UPPER_TRIM_BOUND)), by = .(v_gene, j_gene, v_trim)]
full$productive = 'nonproductive'
full$vj_insert = 0
full$sample_name = 'igor_simulation'

temp_data_filled = merge(temp_data_condensed, full, by = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'productive', 'vj_insert', 'sample_name'), all = TRUE)
temp_data_filled[is.na(count), count := 0]

# merge baseline Igor params
temp_data_filled = merge(temp_data_filled, igor_params, by = c('v_gene', 'j_gene', 'v_trim', 'j_trim'))

# Adjust trimming sites if necessary based on MH ligation
temp_data_adjusted = adjust_trimming_sites_for_ligation_mh(temp_data_filled)

# not going to group genes by common sequence since igor parameters don't do this
temp_data_adjusted$v_gene_group = temp_data_adjusted$v_gene
temp_data_adjusted$j_gene_group = temp_data_adjusted$j_gene

# get gene sequences
whole_nucseq = get_oriented_whole_nucseqs()
temp_data_adjusted = merge(temp_data_adjusted, whole_nucseq[, c('gene', 'v_gene_sequence')], by.x = 'v_gene_group', by.y = 'gene')
temp_data_adjusted = merge(temp_data_adjusted, whole_nucseq[, c('gene', 'j_gene_sequence')], by.x = 'j_gene_group', by.y = 'gene')
setnames(temp_data_adjusted, 'v_gene_sequence', 'v_gene_whole_seq')
setnames(temp_data_adjusted, 'j_gene_sequence', 'j_gene_whole_seq')

# re-condense
cols = c('v_gene_group', 'j_gene_group', 'v_trim', 'j_trim', 'ligation_mh', 'v_trim_prob', 'j_trim_prob', 'v_gene_whole_seq', 'j_gene_whole_seq', 'vj_insert')
temp_data_adjusted_cond = temp_data_adjusted[, sum(count), by = cols]
setnames(temp_data_adjusted_cond, 'V1', 'count')

# get missing observations and sequence context
# Retrieve trim and gene orders based on input types
trims = get_trim_order(TRIM_TYPE)
trim_vars = get_trim_vars(TRIM_TYPE)
genes = get_gene_order(GENE_NAME)

# Filter data based on trim bounds
if (length(trims) > 1){
    temp_data_adjusted_cond = temp_data_adjusted_cond[get(trims[1]) >= LOWER_TRIM_BOUND & get(trims[1]) <= UPPER_TRIM_BOUND & get(trims[2]) >= LOWER_TRIM_BOUND & get(trims[2]) <= UPPER_TRIM_BOUND]
} else {
     temp_data_adjusted_cond = temp_data_adjusted_cond[get(trims[1]) >= LOWER_TRIM_BOUND & get(trims[1]) <= UPPER_TRIM_BOUND]
}

# Extract motifs for observed trims
motif_dataframe = temp_data_adjusted_cond
for (i in seq(length(trims))){
    u_cols = c(paste0(genes[i], '_whole_seq'), paste0(genes[i], '_group'), trims[i])
    subset = unique(motif_dataframe[, ..u_cols])
    subset = subset[, c(paste0(trims[i], '_left_nucs'), paste0(trims[i], '_right_nucs')):= get_nuc_context(get(paste0(genes[i], '_whole_seq')), get(trims[i]))] 
    motif_dataframe = merge(motif_dataframe, subset, by = u_cols)
}

# filter possible sites
# Read only necessary columns and apply filter
frames = fread('https://raw.githubusercontent.com/phbradley/conga/master/conga/tcrdist/db/combo_xcr.tsv', select = c('id', 'region', 'nucseq', 'organism', 'chain', 'frame'))[organism == 'human' & chain == substring(CHAIN_TYPE, 3, 3)]

# Combine operations to reduce redundancy
v = frames[region == 'V', .(v_gene = id, v_frame = frame, v_seq = nucseq)]
j = frames[region == 'J', .(j_gene = id, j_frame = frame, j_seq = nucseq)]

# Merge operations
v$dummy = 1
j$dummy = 1
gene_pairs = merge(v, j, by = 'dummy', allow.cartesian = TRUE)

# Get all trimming sites
trims = data.table(expand.grid(v_trim = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND + 10),
                               j_trim = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND + 10)))
trims$dummy = 1

# Merge genes and trims
all = merge(gene_pairs, trims, by = 'dummy', allow.cartesian = TRUE)[, -c('dummy')]
# Get sequence lengths and subset data to necessary columns
all[, c('v_seq_len', 'j_seq_len') := .(nchar(v_seq), nchar(j_seq))]
cols = c('v_gene', 'j_gene', 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'v_trim', 'j_trim')
all = unique(all[, ..cols])

# Adjust trimming sites based on ligation MH
adjusted_all = adjust_trimming_sites_for_ligation_mh(all)

# Get oriented full sequences and group genes by common features, also subset columns again
cols2 = c('v_gene_group', 'j_gene_group', 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'v_trim', 'j_trim', 'ligation_mh', 'v_gene_sequence', 'j_gene_sequence')
adjusted_all$v_gene_group = adjusted_all$v_gene
adjusted_all$j_gene_group = adjusted_all$j_gene
adjusted_all = merge(adjusted_all, whole_nucseq[, c('gene', 'v_gene_sequence')], by.x = 'v_gene_group', by.y = 'gene')
adjusted_all = merge(adjusted_all, whole_nucseq[, c('gene', 'j_gene_sequence')], by.x = 'j_gene_group', by.y = 'gene')

adjusted_grouped = unique(adjusted_all[, ..cols2])

# Filter by trimming length
adjusted_grouped = adjusted_grouped[v_trim <= UPPER_TRIM_BOUND & j_trim <= UPPER_TRIM_BOUND]
# Get stop positions
adjusted_grouped = get_stop_codon_positions(adjusted_grouped)

# Frame calculations
## j_frame is subtracted because the overall sequence length should have a count of j_frame excess nucleotides
adjusted_grouped[, overall_frame := (v_seq_len + j_seq_len - (v_trim + j_trim) - ligation_mh - j_frame) %% 3]
adjusted_grouped[, frame_type := fifelse(overall_frame == 0, 'In', 'Out')]

# Look for stop codons based on frame
adjusted_grouped = get_stop_positions_with_frame(adjusted_grouped)

# Subset data by columns
cols3 = c('v_gene_group', 'j_gene_group', 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'v_trim', 'j_trim', 'ligation_mh', 'overall_frame', 'frame_type', 'frame_stop')
frame_data = adjusted_grouped[, ..cols3]

# Define frame-related columns
frame_cols = c('frame_type', 'frame_stop', 'overall_frame')

# Remove frame-related columns if present in motif_data
if (any(frame_cols %in% colnames(motif_dataframe))){
    cols = colnames(motif_dataframe)[!(colnames(motif_dataframe) %in% frame_cols)]
    motif_dataframe = motif_dataframe[, ..cols]
}

# Read frame data and filter for possible sites
genes = get_gene_order(GENE_NAME)

# Read frame data and filter for possible sites
possible_sites = frame_data[frame_type == 'Out' | frame_stop == TRUE] 

# Define columns for filtering possible sites
cols = c(paste0(genes, '_group'), 'frame_type', 'overall_frame', 'frame_stop', 'v_trim', 'j_trim', 'ligation_mh')
possible_sites_subset = unique(possible_sites[, ..cols])

# Merge motif data with possible sites subset
cols2 = c(paste0(genes, '_group'), 'v_trim', 'j_trim', 'ligation_mh')
tog = merge(motif_dataframe, possible_sites_subset, by = cols2)

# fill in remaining unobserved, but possible sites 
## Note: we are not including annotations that can be adjusted to a ligation-mh scenario; while these sites are possible, we are ignoring them since we have moved their counts to the ligation-mh scenario

# NOTE: I can't fill in missing possible sites since I don't have IGOR
# parameters for v_trim greater than 14
# motif_data = fill_in_missing_possible_sites(possible_sites_subset, tog, TRIM_TYPE, GENE_NAME)
motif_data = tog

# get scenario weight
motif_data[, total_tcr := sum(count)]
col = c(paste0(genes, '_group'))
motif_data[, paste0(GENE_NAME, '_count') := sum(count), by = col]
motif_data[[paste0('p_', GENE_NAME)]] = motif_data[[paste0(GENE_NAME, '_count')]]/motif_data$total_tcr
motif_data[[paste0('p_', TRIM_TYPE, '_given_', GENE_NAME)]] = motif_data$count/motif_data[[paste0(GENE_NAME, '_count')]]
motif_data$weighted_observation = motif_data[[paste0('p_', TRIM_TYPE, '_given_', GENE_NAME)]]*motif_data[[paste0('p_', GENE_NAME)]]
motif_data$gene_weight_type = 'p_gene_pooled' 


# subset data cols
cols = c('v_gene_group', 'j_gene_group', 'v_trim', 'j_trim', 'ligation_mh', 'v_trim_prob', 'j_trim_prob', 'weighted_observation', 'count', 'total_tcr')
weighted_together_subset = motif_data[, ..cols]

fwrite(weighted_together_subset, filename, sep = '\t')

MODEL_TYPE <<- 'adjusted_igor_norm_ligation-mh'
filename = processed_data_path()
adjusted_norm = weighted_together_subset
adjusted_norm[, j_trim_prob_norm := (j_trim_prob*8)/1]
adjusted_norm[, v_trim_prob_norm := (v_trim_prob*8)/1]
fwrite(adjusted_norm, filename, sep = '\t')

#########################################
## BASELINE IGOR DATA PROCESSING ########
#########################################

MODEL_TYPE <<- 'baseline_igor_ligation-mh'
filename = processed_data_path()

baseline = weighted_together_subset[, -c('v_trim_prob', 'j_trim_prob')]
baseline = merge(baseline, igor_params, by.x = c('v_gene_group', 'j_gene_group', 'v_trim', 'j_trim'), by.y = c('v_gene', 'j_gene', 'v_trim', 'j_trim'))

fwrite(baseline, filename, sep = '\t')

MODEL_TYPE <<- 'baseline_igor_norm_ligation-mh'
filename = processed_data_path()
baseline_norm = baseline
baseline_norm[, j_trim_prob_norm := (j_trim_prob*8)/1]
baseline_norm[, v_trim_prob_norm := (v_trim_prob*8)/1]
fwrite(baseline_norm, filename, sep = '\t')

