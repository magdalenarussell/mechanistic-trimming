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

# LIGATION_PARAMS <<- c(1, 1.25, 1.5, 2)
LIGATION_PARAMS <<- c(1, 1.25, 1.5)

TRIMMING_PROB_MODELS <<- c('motif_two-side-base-count-beyond', 'igor')
MODEL_TYPES <<- c('motif_two-side-base-count-beyond_ligation-mh', 'baseline_igor_norm_ligation-mh')

all_coefs = data.table()

for (param in LIGATION_PARAMS){
    for (trim_model in TRIMMING_PROB_MODELS){
        for (model_type in MODEL_TYPES){
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

            ANNOTATION_TYPE <<- paste0(old_annotation, '_from_', trim_model, '_MHprob', param)

            # Read in model coefficient data 
            coef_path = get_model_coef_file_path(L2)
            coefs = fread(coef_path)

            coefs$model_type = model_type
            coefs$ligation_param = param
            coefs$trimming_prob_model = trim_model

            all_coefs = rbind(coefs, all_coefs)
        }
    }
}

all_coefs[coefficient %like% 'mh' , type:= 'MH parameter']
all_coefs[!(coefficient %like% 'mh'), type:= 'trimming parameter']
subset = all_coefs[, sum(abs(value/log(10))), by = .(type, model_type, ligation_param, trimming_prob_model)]
setnames(subset, 'V1', 'summed_coefficient')

subset_wide = subset %>%
                pivot_wider(names_from = type, values_from = summed_coefficient)%>%
                as.data.table()

subset_wide[model_type %like% 'baseline' & trimming_prob_model %like% 'motif', long_type := paste0('trims drawn from motif + base-count-beyond,\ntraining IGoR + MH model')]
subset_wide[model_type %like% 'motif' & trimming_prob_model %like% 'motif', long_type := paste0('trims drawn from motif + base-count-beyond,\ntraining motif + base-count-beyond + MH model')]
subset_wide[model_type %like% 'baseline' & trimming_prob_model %like% 'ig', long_type := paste0('trims drawn from IGoR,\ntraining IGoR + MH model')]
subset_wide[model_type %like% 'motif' & trimming_prob_model %like% 'ig', long_type := paste0('trims drawn from IGoR,\ntraining motif + base-count-beyond + MH model')]


fwrite(subset_wide, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/scatter_rel_MH_trim_signal/all_coefs.tsv'), sep = '\t')

scatter = ggplot(subset_wide)+
    geom_abline(intercept = 0, slope = 1, size = 3, color = 'gray60', linetype = 'dashed') +
    geom_point(aes(x = get('MH parameter'), y = get('trimming parameter'), color = ligation_param, shape = long_type), size = 10) +
    xlab('magnitude of MH coefficient')+
    xlab('magnitude of trimming coefficients')+
    theme_cowplot(font_family = 'Arial') + 
    labs(color = 'MH simulation parameter', shape = 'Trimming probability model and model trained') +
    background_grid(major = 'xy') + 
    scale_color_distiller(palette = 'Oranges', na.value = 'gray80', direction = +1) + 
    panel_border(color = 'gray60', size = 1.1) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 16), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'), legend.key.width = unit(1.5, 'cm'))  

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/scatter_rel_MH_trim_signal/rel_coefs.pdf')

ggsave(file_name, plot = scatter, width = 15, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
