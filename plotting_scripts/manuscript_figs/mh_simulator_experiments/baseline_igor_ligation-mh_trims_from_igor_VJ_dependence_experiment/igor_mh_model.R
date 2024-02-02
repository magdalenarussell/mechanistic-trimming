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

PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(0)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(0)

MODEL_TYPE <<- 'baseline_igor_norm_ligation-mh'

L2 <<- 'False'

LIGATION_PARAMS <<- c(1, 1.5, 2, 3)
TRIMMING_PROB_MODEL <<- 'igor'

for (param in LIGATION_PARAMS){
    ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
    old_annotation = ANNOTATION_TYPE

    source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))
    source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

    ANNOTATION_TYPE <<- paste0(old_annotation, '_from_', TRIMMING_PROB_MODEL, '_MHprob', param, '_VJ_dependence')
    # Read in model coefficient data 
    coef_path = get_model_coef_file_path(L2)
    coefs = fread(coef_path)

    fwrite(coefs, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/baseline_igor_ligation-mh_trims_from_igor_VJ_dependence_experiment/coefs_MHprob', param, '.tsv'), sep = '\t')

    v_igor_heatmap = plot_igor_coefficient_heatmap_single_group(coefs[coefficient %like% 'v_trim_prob'], with_values = TRUE, limits = c(-0.826, 0.826))+ ggtitle('           V-trimming')

    j_igor_heatmap = plot_igor_coefficient_heatmap_single_group(coefs[coefficient %like% 'j_trim_prob'], with_values = TRUE, limits = c(-0.826, 0.826)) + ggtitle('          J-trimming')

    v_igor_heatmap = v_igor_heatmap + theme(legend.position = 'none') 
    j_igor_heatmap = j_igor_heatmap + theme(legend.position = 'none') 

    # plot mh heatmap
    mh_heatmap = plot_ligation_mh_coefficient_heatmap_single_group(coefs, with_values = TRUE, limits = c(-0.826, 0.826))

    # isolate legend
    legend = get_legend(mh_heatmap) 
    mh_heatmap = mh_heatmap + theme(legend.position = 'none')

    all = align_plots(v_igor_heatmap, j_igor_heatmap, mh_heatmap, legend, align = 'vh', axis = 'lbr')

    first_grid = plot_grid(all[[1]], all[[3]], all[[2]], nrow = 1, rel_widths = c(1, 1, 1), align = 'h')
    together = plot_grid(first_grid, NULL, legend, nrow = 3, rel_heights = c(1, 0.08, 0.2))

    # save plot
    file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/baseline_igor_ligation-mh_trims_from_igor_VJ_dependence_experiment/coef_heatmap_MHprob', param, '.pdf')

    ggsave(file_name, plot = together, width = 18, height = 8, units = 'in', dpi = 750, device = cairo_pdf)
}
