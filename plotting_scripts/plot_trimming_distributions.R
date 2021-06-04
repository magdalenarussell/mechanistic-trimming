source('config/config.R')
args = commandArgs(trailingOnly=TRUE)

TRIM_TYPE <<- args[1]
stopifnot(TRIM_TYPE == 'v_trim')

MOTIF_TYPE <<- 'pnuc_motif' 
stopifnot(MOTIF_TYPE %in% c('pnuc_motif', 'no_pnuc_motif'))

PFM_TYPE <<- 'bounded'
stopifnot(PFM_TYPE %in% c('bounded', 'unbounded'))

NCPU <<- args[2]

REGRESSION_TYPE <<- 'all_subject' 

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 4
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 4

UPPER_TRIM_BOUND <<- 18
LOWER_TRIM_BOUND <<- 2

library(speedglm)
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(NCPU)
library(Biostrings)
library(RhpcBLASctl)
omp_set_num_threads(NCPU)
blas_set_num_threads(NCPU)
library(cowplot)
library(Cairo)

source('scripts/pfm_functions.R')
source('plotting_scripts/plot_dist_functions.R')

actual_trim_data = concatenate_motifs_all_subjects()
predicted_trim_data = fread(get_predicted_dist_file_name())

genes = get_all_genes() 

for (gene_subset in genes){
    plot_predicted_trimming_dists(actual_trim_data, predicted_trim_data, gene_subset)
}

