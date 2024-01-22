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

ANNOTATION_TYPE <<- 'igor_sim_alpha' 

PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(0)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(0)

MODEL_TYPE <<- 'adjusted_igor_norm_ligation-mh'

L2 <<- 'False'

# source(paste0(MOD_PROJECT_PATH, '/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

# Read in model coefficient data 
coef_path = get_model_coef_file_path(L2)
coefs = fread(coef_path)

fwrite(coefs, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/adjusted_igor_norm_ligation-mh_igor_sim_experiment/coefs.tsv'), sep = '\t')

v_igor_heatmap = plot_igor_coefficient_heatmap_single_group(coefs[coefficient %like% 'v_trim_prob'], with_values = FALSE, limits = c(-0.65, 0.65))+ ggtitle('           V-trimming')

j_igor_heatmap = plot_igor_coefficient_heatmap_single_group(coefs[coefficient %like% 'j_trim_prob'], with_values = FALSE, limits = c(-0.65, 0.65)) + ggtitle('          J-trimming')

v_igor_heatmap = v_igor_heatmap + theme(legend.position = 'none') 
j_igor_heatmap = j_igor_heatmap + theme(legend.position = 'none') 

# plot mh heatmap
mh_heatmap = plot_ligation_mh_coefficient_heatmap_single_group(coefs, with_values = FALSE, limits = c(-0.65, 0.65))

# isolate legend
legend = get_legend(mh_heatmap) 
mh_heatmap = mh_heatmap + theme(legend.position = 'none')

all = align_plots(v_igor_heatmap, j_igor_heatmap, mh_heatmap, legend, align = 'vh', axis = 'lbr')

first_grid = plot_grid(all[[1]], all[[3]], all[[2]], nrow = 1, rel_widths = c(1, 1, 1), align = 'h')
together = plot_grid(first_grid, NULL, legend, nrow = 3, rel_heights = c(1, 0.08, 0.2))

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/adjusted_igor_norm_ligation-mh_igor_sim_experiment/coef_heatmap.pdf')

ggsave(file_name, plot = together, width = 18, height = 8, units = 'in', dpi = 750, device = cairo_pdf)
