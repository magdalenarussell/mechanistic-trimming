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
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

MODEL_TYPE <<- 'motif_two-side-base-count-beyond_ligation-mh'

L2 <<- 'True'

source(paste0(MOD_PROJECT_PATH, '/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

filename = processed_data_path()
motif_data = fread(filename)

v_reshaped = data.table(trim_type = 'v_trim', trim_length = motif_data$v_trim, ligation_mh = motif_data$ligation_mh)
j_reshaped = data.table(trim_type = 'j_trim', trim_length = motif_data$j_trim, ligation_mh = motif_data$ligation_mh)
tog_reshaped = data.table(trim_type = 'v_trim + j_trim', trim_length = motif_data$v_trim + motif_data$j_trim, ligation_mh = motif_data$ligation_mh)

reshaped = rbind(v_reshaped, j_reshaped, tog_reshaped)

plot = ggplot(reshaped) +
    geom_point(aes(x = trim_length, y = ligation_mh), size = 4, alpha = 0.5) +
    facet_wrap(~trim_type)+
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 35), axis.text.x=element_text(size = 25), axis.text.y = element_text(size = 25), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    xlab('Trim length') +
    ylab('Ligation MH count')

file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_by_trim_length/trim_length_ligation_mh.pdf')

ggsave(file_name, plot = plot, width = 20, height = 6, units = 'in', dpi = 750, device = cairo_pdf)
