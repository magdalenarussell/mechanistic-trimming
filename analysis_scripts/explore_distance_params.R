source('config/config.R')

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


args = commandArgs(trailingOnly=TRUE)

ANNOTATION_TYPE<<- args[1]
TRIM_TYPE <<- args[2]
trim_types = list.files(path = 'scripts/gene_specific_functions/')
trim_types = str_sub(trim_types, end = -3)
stopifnot(TRIM_TYPE %in% trim_types)

PRODUCTIVITY <<- args[3]

MOTIF_TYPE <<- args[4] 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[5])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- args[6]
GENE_WEIGHT_TYPE <<- args[7]

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[8])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[9])

UPPER_TRIM_BOUND <<- as.numeric(args[10]) 

MODEL_TYPE <<- args[11]
stopifnot(MODEL_TYPE %like% 'distance')

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[12])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')

# Read in model coefficient data 
pwm = get_model_coefficient_data() 

distance_terms = pwm[parameter %like% 'trim_length']

distance_terms$log10_coef = distance_terms$coefficient/log(10)
distance_terms$clean_dist = sapply(distance_terms$parameter, function(x) str_split(x, '_')[[1]][3])
distance_terms$clean_dist = factor(distance_terms$clean_dist, levels = paste0(seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)))

plot = ggplot(distance_terms, aes(x = clean_dist, y = log10_coef)) +
    geom_bar(stat="identity") +
    xlab('Distance') +
    ylab('log10(probability of deletion') +
    theme_cowplot(font_family = 'Arial') + 
    theme(legend.position = "none", text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    ylim(-1.9, 1.9)

path = get_coef_heatmap_file_path()
file_name = paste0(path, '/distance_barplot.pdf')

ggsave(file_name, plot = plot, width = 14, height = 8, units = 'in', dpi = 750, device = cairo_pdf)



