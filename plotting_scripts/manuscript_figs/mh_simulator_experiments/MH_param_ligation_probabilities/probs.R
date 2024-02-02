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
# LIGATION_PARAMS <<- c(1, 1.25, 1.5, 2.0, 3.0)
LIGATION_PARAMS <<- c(1, 1.25, 1.5)


source(paste0(MOD_PROJECT_PATH,'/analysis_scripts/ligation-mh_signal_simulator_scripts/ligation-mh_simulator_functions.R'))

#TODO remove this when doing new analysis
get_ligation_probabilities <- function(lig_param, possible_ligs){
    stopifnot(lig_param >= 1)
    starting_prob = 0.25
    lig_probs = rep(starting_prob, length(possible_ligs))
    names(lig_probs) = possible_ligs    
    if (lig_param > 1){
        scale = lig_param * possible_ligs[possible_ligs > 0] 
    } else {
        scale = lig_param
    }
    lig_probs[names(lig_probs) > 0] = lig_probs[names(lig_probs) > 0]*scale
    lig_probs[names(lig_probs) == 0] = lig_probs[names(lig_probs) == 0]* (1/lig_param)
    lig_probs['no_ligation'] = 0.25 - 0.1
    prob_sum = sum(lig_probs)
    lig_probs = lig_probs/prob_sum
    return(lig_probs)
}

lig_probs = data.table()

for (lig in LIGATION_PARAMS){
    temp = get_ligation_probabilities(lig, seq(0, 8))
    temp_df = data.table(ligation_mh = names(temp), prob = temp)
    temp_df$ligation_param = lig
    lig_probs = rbind(lig_probs, temp_df)
}

lig_probs$ligation_mh = factor(lig_probs$ligation_mh,levels = c('no_ligation', seq(0, 8)))
# lig_probs$ligation_param = factor(lig_probs$ligation_param, levels = c(3.0, 2.0, 1.5, 1.25, 1.0))

scatter = ggplot(lig_probs)+
    geom_point(aes(x = ligation_mh, y = log10(prob), color = ligation_param), size = 10) +
    geom_line(aes(x = ligation_mh, y = log10(prob), group = ligation_param, color = ligation_param), size = 5) +
    xlab('Amount of MH involved in ligation') +
    ylab('log10(Ligation probability)') +
    theme_cowplot(font_family = 'Arial') + 
    labs(color = 'MH simulation parameter') +
    background_grid(major = 'xy') + 
    scale_color_distiller(palette = 'Oranges', na.value = 'gray80', direction = +1) + 
    panel_border(color = 'gray60', size = 1.1) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 16), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'), legend.key.width = unit(1.5, 'cm'))  

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/MH_param_ligation_probabilities/lig_params_log.pdf')

ggsave(file_name, plot = scatter, width = 15, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

scatter = ggplot(lig_probs)+
    geom_point(aes(x = ligation_mh, y = prob, color = ligation_param), size = 10) +
    geom_line(aes(x = ligation_mh, y = prob, group = ligation_param, color = ligation_param), size = 5) +
    xlab('Amount of MH involved in ligation') +
    ylab('Ligation probability') +
    theme_cowplot(font_family = 'Arial') + 
    labs(color = 'MH simulation parameter') +
    background_grid(major = 'xy') + 
    scale_color_distiller(palette = 'Oranges', na.value = 'gray80', direction = +1) + 
    panel_border(color = 'gray60', size = 1.1) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 16), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'), legend.key.width = unit(1.5, 'cm'))  

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/MH_param_ligation_probabilities/lig_params.pdf')

ggsave(file_name, plot = scatter, width = 15, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

