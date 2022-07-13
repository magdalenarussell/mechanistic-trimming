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


args = commandArgs(trailingOnly=TRUE)

ANNOTATION_TYPE <<- args[1]

TRIM_TYPE <<- 'v_trim' 

PRODUCTIVITY <<- args[2]

MOTIF_TYPE <<- args[3] 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[4])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[5]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 0 
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 0

UPPER_TRIM_BOUND <<- as.numeric(args[6]) 
LOWER_TRIM_BOUND <<- 2 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA

MODEL_TYPE <<- 'distance'

TYPE <<- args[7]
stopifnot(TYPE %like% 'v_gene_family_loss')

source('scripts/data_compilation_functions.R')
source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/model_evaluation_functions.R')

if (TYPE == 'v_gene_family_loss'){
    cluster_count = 3
    combine_by_terminal = TRUE
    full_sequence = FALSE
    align = FALSE
} else if (TYPE == 'full_v_gene_family_loss'){
    cluster_count = 4 
    combine_by_terminal = FALSE
    full_sequence = TRUE
    align = TRUE 
}

v_families = get_gene_families(cluster_count, combine_by_terminal, full_sequence, align)$cluster_data

TRIM_TYPE <<- 'j_trim'
GENE_NAME <<- 'j_gene'

source('scripts/data_compilation_functions.R')
source('scripts/model_evaluation_functions.R')

j_families = get_gene_families(cluster_count = 1, combine_by_terminal, full_sequence, align)$cluster_data

all_seqs = rbind(unique(v_families[, c('gene', 'terminal_seq', 'clusters_grouped')]), unique(j_families[, c('gene', 'terminal_seq', 'clusters_grouped')]), fill = TRUE)

seq_list = DNAStringSet(all_seqs$terminal_seq)

require(ape)
require(DECIPHER)

if (isTRUE(align)) {
    seq_list = AlignSeqs(seq_list)
}

dists = DistanceMatrix(seq_list)

colnames(dists) = all_seqs$gene
row.names(dists) = all_seqs$gene
dist_format = as.dist(dists) 
   
clusters = hclust(dist_format)
# clusters_grouped = cutree(clusters, k = 1)
all_seqs[substring(gene, 4, 4) == 'J', clusters_grouped := cluster_count +1]
clusters_grouped = all_seqs$clusters_grouped
names(clusters_grouped) = all_seqs$gene

    
require(RColorBrewer)
colors = brewer.pal(cluster_count +1, 'Set2')

TRIM_TYPE <<- 'v_trim'
path = get_analysis_path()
filename2 = file.path(path, paste0(TYPE, '_unrooted_tree_v_j_comparison.pdf'))

cairo_pdf(filename2, width = 12, height = 12, fallback_resolution = 750)
plot(as.phylo(clusters), type = 'unrooted', tip.color = colors[clusters_grouped], no.margin = TRUE)
dev.off()
