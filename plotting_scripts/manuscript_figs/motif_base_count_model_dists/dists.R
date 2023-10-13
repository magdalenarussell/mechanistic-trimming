source('mechanistic-trimming/config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(ggplot2)
library(cowplot)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

ANNOTATION_TYPE <<- 'igor_alpha' 

PARAM_GROUP <<- 'nonproductive_v-j_trim'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

MODEL_TYPE <<- 'motif_two-side-base-count-beyond'

L2 <<- 'True'

source(paste0(MOD_PROJECT_PATH, '/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

# Read in dist data
predicted_trims_path = get_model_predictions_file_path(L2)
predicted_trims = fread(predicted_trims_path) 
predicted_trims[, gene_pair := paste0(v_gene_group, ' + ', j_gene_group)]

# subset to top genes
predicted_trims[, p_gene_pair := sum(count)/total_tcr, by = .(gene_pair)]
subset = predicted_trims[p_gene_pair > 3.05e-3]

# Create a trimming distribution plot for each gene
v_plot = plot_predicted_paired_trimming_dists_single_group(subset, row_var = 'j_trim', color = '#01665e') 
# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/motif_base_count_model_dists/v_dist_gallery.pdf')
ggsave(file_name, plot = v_plot, width = 30, height = 35, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)

j_plot = plot_predicted_paired_trimming_dists_single_group(subset, row_var = 'v_trim', color = '#01665e') 
# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/motif_base_count_model_dists/j_dist_gallery.pdf')
ggsave(file_name, plot = j_plot, width = 30, height = 35, units = 'in', dpi = 750, device = cairo_pdf, limitsize = FALSE)


