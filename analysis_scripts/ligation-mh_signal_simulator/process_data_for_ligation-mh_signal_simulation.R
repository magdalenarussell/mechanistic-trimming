source('mechanistic-trimming/config/config.R')
source('mechanistic-trimming/config/file_paths.R')

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

args = commandArgs(trailingOnly=TRUE)

ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
PARAM_GROUP <<- args[1]
stopifnot(PARAM_GROUP %in% c('nonproductive_v-j_trim_ligation-mh', 'both_v-j_trim_ligation-mh'))
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))
NCPU <<- as.numeric(args[2])
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond_ligation-mh'

TRIMMING_PROB_TYPE <<- args[3]
stopifnot(TRIMMING_PROB_TYPE %in% c('igor', 'motif_two-side-base-count-beyond', 'uniform', 'mh_adjusted_motif-two-side-base-count-beyond'))

LIGATION_MH_PARAM <<- as.numeric(args[4])

source(paste0(MOD_PROJECT_PATH,'/analysis_scripts/ligation-mh_signal_simulator/ligation-mh_simulator_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# total number of sequences
total = 2500000

# get V and J gene choice probs
gene_choice_probs = get_igor_gene_usage_params()
whole_nucseq = get_oriented_whole_nucseqs()
gene_choice_probs$v_choice = merge(gene_choice_probs$v_choice, whole_nucseq[, c('gene', 'v_gene_sequence')], by.x = 'v_gene', by.y = 'gene') 
gene_choice_probs$j_choice = merge(gene_choice_probs$j_choice, whole_nucseq[, c('gene', 'j_gene_sequence')], by.x = 'j_gene', by.y = 'gene') 

# get Vtrim and Jtrim probs
trim_probs = get_trimming_probs(TRIMMING_PROB_TYPE)

# get all ligation probabilities
all_lig_probs = get_ligation_probabilities(LIGATION_MH_PARAM, seq(0, 15))
all_lig_probs = all_lig_probs[names(all_lig_probs) != 'no_ligation']

# get all possible configurations
configs = read_frames_data()

set.seed(123) 

registerDoParallel(cores=NCPU)
sim_data = foreach(subset = seq(total/100), .combine=rbind) %dopar% {
    subset_sim = data.table()
    for (i in seq(100)){
        # print(paste0('starting sim ', i))
        # sample a V and J gene
        temp = cbind(gene_choice_probs$v_choice[sample(.N, 1, prob = v_gene_prob)],
                     gene_choice_probs$j_choice[sample(.N, 1, prob = j_gene_prob)]) 

        # sample a V and J gene trimming amount
        v_possible_trims = trim_probs$v_trim[v_gene == temp$v_gene]
        j_possible_trims = trim_probs$j_trim[j_gene == temp$j_gene]

        if (nrow(v_possible_trims) == 0 | nrow(j_possible_trims) == 0){
            next
        }

        if (all(unique(v_possible_trims$v_trim_prob) == 0) | all(unique(j_possible_trims$j_trim_prob) == 0)){
            next
        }

        cols = c('v_gene', 'v_gene_prob', 'v_gene_sequence', 'j_gene', 'j_gene_prob', 'j_gene_sequence')
        temp = temp[, ..cols]

        ## Create a probability vector where lower values have higher probability
        ## For example, using linearly decreasing weights
        temp_vt = v_possible_trims[sample(.N, 1, prob=v_trim_prob)]
        temp_jt = j_possible_trims[sample(.N, 1, prob=j_trim_prob)]

        temp = merge(temp, temp_vt, by = 'v_gene')
        temp = merge(temp, temp_jt, by = 'j_gene')

        temp_configs = configs[v_gene == temp$v_gene & j_gene == temp$j_gene & v_trim == temp$v_trim & j_trim == temp$j_trim]

        lig_probs = all_lig_probs[names(all_lig_probs) %in% unique(temp_configs$ligation_mh)]

        lig_draw = sample(names(lig_probs), 1, prob = lig_probs)

        temp$ligation_mh = as.numeric(lig_draw)
        subset_sim = rbind(subset_sim, temp)
    }
    subset_sim
}
stopImplicitCluster()

cols = colnames(sim_data)[colnames(sim_data) != 'ligation_attempt']

condensed_sim = sim_data[, .N, by = cols]
setnames(condensed_sim, 'N', 'count')

condensed_sim$vj_insert = 0

motif_data = get_all_nuc_contexts(condensed_sim, 'mh_sim', gene_type = GENE_NAME, trim_type = TRIM_TYPE)

# Because there is no selection effect in these simulated data, I am not
# restricting to nonproductive sequences! 
# fill in missing sites
cols2 = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'ligation_mh')
tog = merge(motif_data, configs, by = cols2)

if (grepl('both', PARAM_GROUP)){
    filled_motif_data = fill_in_missing_possible_sites(configs, tog, trim_type = TRIM_TYPE, gene_type = GENE_NAME)
} else if (grepl('nonprod', PARAM_GROUP)){
    filled_motif_data = filter_motif_data_for_possible_sites(tog)
}

processed = inner_aggregation_processing(filled_motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE)
processed = subset_processed_data(processed, trim_type = TRIM_TYPE, gene_type = GENE_NAME)

old_ANNOTATION_TYPE <<- ANNOTATION_TYPE
ANNOTATION_TYPE <<- paste0(old_ANNOTATION_TYPE, '_from_', TRIMMING_PROB_TYPE, '_MHprob', LIGATION_MH_PARAM)

filename = processed_data_path()
fwrite(processed, filename, sep = '\t')

MODEL_TYPE <<- 'baseline_igor_norm_ligation-mh'
LEFT_NUC_MOTIF_COUNT <<- 0
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 0 

if (length(colnames(processed)[colnames(processed) %like% 'ligation']) > 1){
    processed = processed[, -23]
}

baseline = merge(processed, trim_probs$v_trim, by = c('v_gene', 'v_trim'))
baseline = merge(baseline, trim_probs$j_trim, by = c('j_gene', 'j_trim'))

filename = processed_data_path()

baseline[, j_trim_prob_norm := (j_trim_prob*8)/1]
baseline[, v_trim_prob_norm := (v_trim_prob*8)/1]
fwrite(baseline, filename, sep = '\t')
