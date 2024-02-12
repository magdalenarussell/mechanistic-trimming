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

LIGATION_PARAMS <<- c(0, 1, 10)
TRIMMING_PROB_MODEL <<- c('igor', 'motif_two-side-base-count-beyond')
ANNOTATION_SUBTYPE <<- c('always-ligate')

all_data = data.table()

for (trim_model in TRIMMING_PROB_MODEL){
    for (subtype in ANNOTATION_SUBTYPE){
        for (param in LIGATION_PARAMS){
            ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
            old_annotation = ANNOTATION_TYPE

            source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))
            source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

            ANNOTATION_TYPE <<- paste0(old_annotation, '_from_', trim_model, '_MHprob', param, '_', subtype)

            # read lrt data
            lrt_path = get_LRT_file_path(L2)
            if (!file.exists(lrt_path)){
                next
            }
            lrts = fread(lrt_path)
            lrts$MH_param = param
            lrts$trimming_prob_model = trim_model

            all_data = rbind(all_data, lrts)
        }
    }
}

fwrite(all_data, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/lrt/lrt_data.tsv'), sep = '\t')

all_data[, MH_param_factor := factor(MH_param)]
all_data[MH_param == 0, MH_param := 0.01]
all_data[trimming_prob_model == 'igor', trim_prob_long := paste0('Simulator trims drawn from\nIGoR model')]
all_data[trimming_prob_model != 'igor', trim_prob_long := paste0('Simulator trims drawn from\nmotif + base-count model')]

insig_thresh = qchisq(0.05, df = 1, lower.tail = FALSE)

scatter = ggplot(all_data)+
    geom_hline(yintercept = insig_thresh, size = 3, color = 'gray60', linetype = 'dashed') +
    geom_point(aes(y = LRT_test_statistic, x = MH_param_factor, color = MH_param), size = 10) +
    facet_wrap(~trim_prob_long, ncol = 2) +
    xlab('MH simulation parameter')+
    ylab('LRT test statistic') +
    theme_cowplot(font_family = 'Arial') + 
    labs(color = 'log(LRT p-value)') +
    background_grid(major = 'xy') + 
    # scale_color_distiller(palette = 'Oranges', direction = +1) +
    scale_color_distiller(palette = 'Oranges', direction = +1, trans = "log10", oob = scales::oob_squish, limits = c(1e-2, 10)) +
    panel_border(color = 'gray60', size = 1.1) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 16), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'), legend.key.width = unit(1.5, 'cm'), legend.position = 'none')  

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/lrt/rel_logLik.pdf')

ggsave(file_name, plot = scatter, width = 15, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
