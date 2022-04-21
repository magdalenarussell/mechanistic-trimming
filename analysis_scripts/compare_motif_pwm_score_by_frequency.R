source('config/config.R')

library(ggplot2)
library(cowplot)
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)


args = commandArgs(trailingOnly=TRUE)

ANNOTATION_TYPE <<- args[1]
stopifnot(ANNOTATION_TYPE %in% c('igor', 'parsimony'))

TRIM_TYPE <<- args[2]
stopifnot(TRIM_TYPE == 'v_trim')

PRODUCTIVITY <<- args[3]

MOTIF_TYPE <<- args[4] 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[5])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- args[6]

GENE_WEIGHT_TYPE <<- args[7]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[8])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[9])

UPPER_TRIM_BOUND <<- as.numeric(args[10]) 

MODEL_TYPE <<- args[11]
stopifnot(MODEL_TYPE %like% 'motif')

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[12])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')
source('analysis_scripts/pwm_profile_functions.R')

# Compile data for all subjects
motif_data = aggregate_all_subject_data()
positions = get_positions()
cols = c('gene', 'trim_length', 'motif', positions)
condensed = unique(motif_data[,..cols]) 

# Read in model coefficient data 
pwm = get_model_coefficient_data() 

# calculate pwm score by gene and trim length
condensed[, pwm_score := unlist(lapply(motif, function(x) as.numeric(get_pwm_score(pwm, x, positions))))]

# merge
together = merge(motif_data, condensed, by = colnames(condensed)[!colnames(condensed) == 'pwm_score'])

# get trim frequencies
freqs = together[, sum(weighted_observation), by = .(motif, pwm_score)]
setnames(freqs, 'V1', 'empirical_trim_prob')

# get trim frequencies
obs_freqs = condensed[, .N, by = .(motif, pwm_score)]
obs_freqs[, N := N/sum(N)]
setnames(obs_freqs, 'N', 'empirical_frequency')

freqs_tog = merge(obs_freqs, freqs)

plot = ggplot(freqs_tog, aes(x = log10(empirical_trim_prob), y = pwm_score/log(10))) +
    geom_point(size = 4, alpha = 0.5) +
    geom_smooth(method = 'lm', size = 2) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('log10(empirical motif trimming frequency)') +
    ylab(paste0('log10(PWM score)\nfrom ', MODEL_TYPE)) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank()) +
    background_grid(major = 'xy') +
    panel_border(color = 'gray60', size = 1.5) 

path = get_pwm_profile_plot_file_path(0)    
ggsave(paste0(path, '/pwm_score_by_motif_empirical_trim_freq.pdf'), plot, width = 8, height = 8, units = 'in', dpi = 750, device = cairo_pdf)

plot = ggplot(freqs_tog, aes(x = log10(empirical_frequency), y = pwm_score/log(10))) +
    geom_point(size = 4, alpha = 0.5) +
    geom_smooth(method = 'lm', size = 2) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('log10(empirical motif frequency)') +
    ylab(paste0('log10(PWM score)\nfrom ', MODEL_TYPE)) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank()) +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

path = get_pwm_profile_plot_file_path(0)    
ggsave(paste0(path, '/pwm_score_by_motif_empirical_freq.pdf'), plot, width = 8, height = 8, units = 'in', dpi = 750, device = cairo_pdf)





