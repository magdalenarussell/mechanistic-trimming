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

ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- 2
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_ligation-mh'

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# total number of sequences
total = 1600000

# get V and J gene sequences
whole_nucseqs = get_whole_nucseqs()
vgenes = whole_nucseqs[gene %like% 'TRAV']
colnames(vgenes) = c('v_gene', 'v_gene_sequence')
jgenes = whole_nucseqs[gene %like% 'TRAJ']
colnames(jgenes) = c('j_gene', 'j_gene_sequence')

# simulate single sequence
set.seed(123) 

sim_data = data.table()

for (i in seq(total)){
    # sample a V and J gene
    temp = cbind(vgenes[sample(.N, 1)],
                 jgenes[sample(.N, 1)])

    # sample a V and J gene trimming amount
    v_possible_trims = seq(-2, 16)
    j_possible_trims = seq(-2, 16)

    try = 0
    while(!('ligation_mh' %in% colnames(temp))) {
        try = try + 1
        print(paste0(try, ' sampling trims'))

        ## Create a probability vector where lower values have higher probability
        ## For example, using linearly decreasing weights
        v_weights = rev(seq_along(v_possible_trims))
        v_weights = v_weights^4 / sum(v_weights^4)
        if (length(v_weights) == 1){
            vt = unique(v_possible_trims)
        } else {
            vt = sample(v_possible_trims, 1, prob = v_weights)
        }

        j_weights = rev(seq_along(j_possible_trims))
        j_weights = j_weights^4 / sum(j_weights^4)
        if (length(j_weights) == 1){
            jt = unique(j_possible_trims)
        } else {
            jt = sample(j_possible_trims, 1, prob = j_weights)
        }

        ## Filter the vector to remove numbers less than or equal to the sampled number
        v_possible_trims = v_possible_trims[v_possible_trims >= vt]
        j_possible_trims = j_possible_trims[j_possible_trims >= jt]

        temp$v_trim = vt
        temp$j_trim = jt

        # get possible ligation MH
        overlaps = seq(1, 15)
        possible_ligs = c()
        for (o in overlaps){
            lig_mh_count = get_possible_ligation_mh_fixed_trim(temp, overlap_count=o)
            possible_ligs = unique(c(possible_ligs, lig_mh_count))
        }

        ordered_possible_ligs = possible_ligs[order(-possible_ligs)]

        lig_probs = rep(0.35, length(possible_ligs))
        lig_probs = lig_probs ^ (1/(ordered_possible_ligs+1))

        for (m in seq(length(lig_probs))){
            print(paste0(try, ' sampling lig for ', ordered_possible_ligs[m]))
            p = lig_probs[m]
            lig_draw = sample(c(TRUE, FALSE), 1, prob = c(p, 1-p))
            if (isTRUE(lig_draw)){
                temp$ligation_mh = ordered_possible_ligs[m]
            }
        }
    }
    
    sim_data = rbind(sim_data, temp)
}


temp[, c("v_gene.j_trimmed", "j_trimmed.v_gene", "j_gene.v_trimmed", "v_trimmed.j_gene") := get_overlapping_regions_ligation_mh(v_gene_sequence, j_gene_sequence, v_trim, j_trim)]

# get MH
temp = get_mh_dataframe(temp, aligning_trim = 'j_trim', aligning_gene = 'v_gene')
temp = get_mh_dataframe(temp, aligning_trim = 'v_trim', aligning_gene = 'j_gene')

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
