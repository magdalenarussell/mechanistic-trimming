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

LIGATION_PARAMS <<- c(0, 0.1, 1, 10)
TRIMMING_PROB_MODEL <<- c('igor', 'motif_two-side-base-count-beyond')
ANNOTATION_SUBTYPE <<- c('always-ligate', 'always-ligate_uniformVJ', 'no-sequential-trim', 'no-sequential-trim_uniformVJ')

for (trim_model in TRIMMING_PROB_MODEL){
    for (subtype in ANNOTATION_SUBTYPE){
        for (param in LIGATION_PARAMS){
            ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
            old_annotation = ANNOTATION_TYPE

            source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))
            source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

            ANNOTATION_TYPE <<- paste0(old_annotation, '_from_', trim_model, '_MHprob', param, '_', subtype)
            # Read in model coefficient data 
            coef_path = get_model_coef_file_path(L2)
            if (!file.exists(coef_path)){
                assign(paste0('first_grid_', param), NULL)

                next
            }
            coefs = fread(coef_path)

            fwrite(coefs, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/baseline_igor_ligation-mh_trims_always-ligate_no-sequential_experiment/coefs_MHprob', param, '_', subtype, '_', trim_model, '_trim-model.tsv'), sep = '\t')

            v_igor_heatmap = plot_igor_coefficient_heatmap_single_group(coefs[coefficient %like% 'v_trim_prob'], with_values = TRUE, limits = c(-0.832, 0.832))+ ggtitle('           V-trimming')

            j_igor_heatmap = plot_igor_coefficient_heatmap_single_group(coefs[coefficient %like% 'j_trim_prob'], with_values = TRUE, limits = c(-0.832, 0.832)) + ggtitle('          J-trimming')

            v_igor_heatmap = v_igor_heatmap + theme(legend.position = 'none') 
            j_igor_heatmap = j_igor_heatmap + theme(legend.position = 'none') 

            # plot mh heatmap
            mh_heatmap = plot_ligation_mh_coefficient_heatmap_single_group(coefs, with_values = TRUE, limits = c(-0.832, 0.832))

            # isolate legend
            legend = get_legend(mh_heatmap) 
            mh_heatmap = mh_heatmap + theme(legend.position = 'none')

            all = align_plots(v_igor_heatmap, j_igor_heatmap, mh_heatmap, legend, align = 'vh', axis = 'lbr')

            assign(paste0('first_grid_', param), plot_grid(all[[1]], all[[3]], all[[2]], nrow = 1, rel_widths = c(1, 1, 1), align = 'h'))
        }

        names = paste0('ligation-MH prob effect = ', LIGATION_PARAMS)
        # Create a new vector of length 2 * length(s) - 1
        names2 <- vector("character", length = 2 * length(names))

        # Fill in the original elements
        names2[seq(1, by = 2, length.out = length(names))] = names 

        # Fill in NA or "" at the alternate positions
        names2[seq(2, by = 2, length.out = length(names))] = ""

        together = plot_grid(NULL, first_grid_0, NULL, first_grid_0.1, NULL, first_grid_1, NULL, first_grid_10, NULL, legend, nrow = 10, rel_heights = c(0.10, 1, 0.10, 1, 0.10, 1, 0.10, 1, 0.10, 0.2), labels = c(names2, "", ""), label_size = 30, hjust = 0, label_x = 0.01)

        # save plot
        file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/baseline_igor_ligation-mh_trims_always-ligate_no-sequential_experiment/coef_heatmap_', subtype, '_', trim_model, '_trim-model.pdf')

        ggsave(file_name, plot = together, width = 18, height = 26, units = 'in', dpi = 750, device = cairo_pdf)
    }
}
