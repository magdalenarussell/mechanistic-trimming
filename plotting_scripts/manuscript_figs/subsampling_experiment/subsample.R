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

MODEL_TYPE <<- 'motif_two-side-base-count-beyond_interior-mh-count'

L2 <<- 'True'

source(paste0(MOD_PROJECT_PATH, '/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

path = get_subsample_file_path(L2)
files = list.files(path, full.names = TRUE, pattern = '*tsv') 
results = data.table()
for (file in files){
    temp = fread(file)
    basename = strsplit(file, "/")[[1]][length(strsplit(file, "/")[[1]])]
    basename = strsplit(basename, 'data_prop')[[1]][2]
    basename = strsplit(basename, '_coefs_')[[1]][1]
    temp$prop = as.numeric(basename)
    results = rbind(results, temp)
}

results$log_coef = results$value/log(10)
results = results[!is.na(log_coef)]
results[, q25 := quantile(log_coef, 0.25), by = .(coefficient, base, position, side, trim_type, prop)]
results[, q75 := quantile(log_coef, 0.75), by = .(coefficient, base, position, side, trim_type, prop)]
results[, median := quantile(log_coef, 0.5), by = .(coefficient, base, position, side, trim_type, prop)]
results[, variance := var(log_coef), by = .(coefficient, base, position, side, trim_type, prop)]

cols = c('coefficient', 'base', 'position', 'side', 'trim_type', 'prop', 'q25', 'q75', 'median', 'variance')
condensed = unique(results[, ..cols])

condensed$prop = factor(condensed$prop, levels = rev(unique(condensed$prop)))

condensed[coefficient == 'mh_count', nice_param := paste0(side, ' MH count, ', substring(position,8,8), ' nt overlap')]
condensed[coefficient == 'base_count', nice_param := paste0(side, ' ', base, ' base-count beyond, ', toupper(substring(trim_type, 1, 1)), '-gene')]
condensed[coefficient == 'motif', nice_param := paste0(side, ' motif pos ', substring(position, 4, 4), ', ', base, ' nt, ', toupper(substring(trim_type, 1, 1)), '-gene')]

condensed$nice_param = factor(condensed$nice_param, levels = condensed[prop == 0.5][order(variance)]$nice_param)

label_data = condensed[prop == 0.05]
color_palette = set_color_palette(unique(condensed$coefficient))

require(ggrepel)
plot = ggplot(condensed[order(match(nice_param, levels(condensed$nice_param)))]) +
    geom_line(aes(x = prop, y = variance, color = coefficient, group = nice_param), size = 4, alpha = 0.8) +
    geom_point(aes(x = prop, y = variance, group = nice_param, color = coefficient), size = 8) +
    geom_text_repel(data = label_data, aes(y = variance, x = prop, label = nice_param, color = coefficient), nudge_x = 0.2, fontface = "bold", size = 11, direction = 'y', hjust = 0, point.padding = 1, max.overlaps = Inf, lineheight = 0.8) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('\nProportion of full training data set') +
    ylab('Variance of log10(probability of deletion)\n') +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    theme(legend.position = 'none', text = element_text(size = 45), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 40), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
    scale_x_discrete(expand = expansion(add = c(0.1, 4))) +
    scale_color_manual(values = color_palette)


file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/subsampling_experiment/subsample.pdf')

ggsave(file_name, plot = plot, width = 25, height = 22, units = 'in', dpi = 750, device = cairo_pdf)
