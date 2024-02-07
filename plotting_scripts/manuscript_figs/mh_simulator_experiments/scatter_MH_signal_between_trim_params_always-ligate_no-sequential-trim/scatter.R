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

L2 <<- 'False'

# LIGATION_PARAMS <<- c(1, 1.5, 2, 3)
LIGATION_PARAMS <<- c(0, 0.1, 1, 10)


TRIMMING_PROB_MODELS <<- c('motif_two-side-base-count-beyond', 'igor')

MODEL_TYPES <<- c('motif_two-side-base-count-beyond_ligation-mh', 'baseline_igor_norm_ligation-mh')
# MODEL_TYPES <<- c('baseline_igor_norm_ligation-mh')

subtype = c('always-ligate', 'no-sequential-trim')
all_coefs = data.table()

for (param in LIGATION_PARAMS){
    for (trim_model in TRIMMING_PROB_MODELS){
        for (model_type in MODEL_TYPES){
            for (sub in subtype){
                if (model_type %like% 'motif'){
                    # 5' motif nucleotide count
                    LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
                    # 3' motif nucleotide count
                    RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)
                } else {
                    # 5' motif nucleotide count
                    LEFT_NUC_MOTIF_COUNT <<- as.numeric(0)
                    # 3' motif nucleotide count
                    RIGHT_NUC_MOTIF_COUNT <<- as.numeric(0)
                }

                ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'

                old_annotation = ANNOTATION_TYPE

                MODEL_TYPE <<- model_type

                source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))
                source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

                ANNOTATION_TYPE <<- paste0(old_annotation, '_from_', trim_model, '_MHprob', param, '_', sub)

                # Read in model coefficient data 
                coef_path = get_model_coef_file_path(L2)
                if (!file.exists(coef_path)){
                    next
                }
 
                coefs = fread(coef_path)

                coefs$model_type = model_type
                coefs$ligation_param = param
                coefs$trimming_prob_model = trim_model
                coefs$simulation_type = sub
                all_coefs = rbind(coefs, all_coefs)
            }
        }
    }
}

subset = all_coefs[coefficient %like% 'mh']
subset$value = subset$value/log(10)

cols = c('coefficient', 'ligation_param', 'trimming_prob_model', 'model_type', 'value', 'simulation_type')
subset_wide = subset[, ..cols] %>%
                pivot_wider(names_from = trimming_prob_model, values_from = value)%>%
                as.data.table()
subset_wide = unique(subset_wide)
fwrite(subset_wide, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/scatter_MH_signal_between_trim_params_always-ligate_no-sequential-trim/all_coefs.tsv'), sep = '\t')

subset_wide$shape_var = paste0(subset_wide$simulation_type, ', ', subset_wide$model_type)
subset_wide$shape_var = str_replace(subset_wide$shape_var, 'baseline_igor_norm_ligation-mh', 'baseline-igor + ligation-mh')
subset_wide$shape_var = str_replace(subset_wide$shape_var, 'motif_two-side-base-count-beyond_ligation-mh', 'motif + base-count + ligation-mh')

# Adjust zero values to a small, non-zero value
subset_wide[ligation_param == 0, ligation_param := 0.01] 

scatter = ggplot(subset_wide)+
    geom_abline(intercept = 0, slope = 1, size = 3, color = 'gray60', linetype = 'dashed') +
    geom_point(aes(x = get('motif_two-side-base-count-beyond'), y = get('igor'), color = ligation_param, shape = shape_var), size = 10) +
    xlab('MH coefficient for trims sampled according to\nmotif + base-count probabilities')+
    ylab('MH coefficient for trims sampled according to\nIGoR probabilities')+
    theme_cowplot(font_family = 'Arial') + 
    labs(color = 'MH simulation parameter', shape = 'Simulation type, Model trained') +
    background_grid(major = 'xy') + 
    scale_color_distiller(palette = 'Oranges', direction = +1, trans = "log10", oob = scales::oob_squish, limits = c(1e-2, 10)) + 
    panel_border(color = 'gray60', size = 1.1) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 16), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'), legend.key.width = unit(1.5, 'cm'))  

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/scatter_MH_signal_between_trim_params_always-ligate_no-sequential-trim/rel_MH_coefs.pdf')

ggsave(file_name, plot = scatter, width = 18, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
