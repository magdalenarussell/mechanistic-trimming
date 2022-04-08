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

TRIM_TYPE <<- args[2] 

PRODUCTIVITY <<- args[3]

MOTIF_TYPE <<- args[4] 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[5])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[6]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 0 
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 0

UPPER_TRIM_BOUND <<- as.numeric(args[7]) 
LOWER_TRIM_BOUND <<- 2 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA

MODEL_TYPE <<- 'distance'

TYPE <<- args[8]
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
    cluster_count = 5 
    combine_by_terminal = FALSE
    full_sequence = TRUE
    align = TRUE 
}

v_families = get_gene_families(cluster_count, combine_by_terminal, full_sequence, align)

clusters_grouped = v_families$cluster_data$clusters_grouped
names(clusters_grouped) = v_families$cluster_data$gene

    
require(RColorBrewer)
colors = brewer.pal(cluster_count, 'Set2')

path = get_analysis_path()
filename2 = file.path(path, paste0(TYPE, '_unrooted_tree.pdf'))

cairo_pdf(filename2, width = 12, height = 12, fallback_resolution = 750)
plot(v_families$tree, type = 'unrooted', tip.color = colors[clusters_grouped], no.margin = TRUE)
dev.off()
