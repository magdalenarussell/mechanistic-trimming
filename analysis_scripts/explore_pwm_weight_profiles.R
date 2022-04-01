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

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)){
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

wide = condensed[, c('gene', 'trim_length', 'pwm_score')] %>% 
    pivot_wider(names_from = 'trim_length', values_from = 'pwm_score') %>%
    as.data.table()

cluster_plot = determine_cluster_count_plot(wide)
cluster_data = cluster_data(wide, 5, 'pwm_score')

cluster_data$trim_length = as.numeric(cluster_data$trim_length)
cluster_data = cluster_data[order(trim_length)]
lines = ggplot(cluster_data) +
    geom_line(aes(x = trim_length, y = pwm_score/log(10), group = gene, color = as.factor(k_means_cluster)), size = 3, alpha = 0.5) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Trim length') +
    ylab('-log10(PWM score)') +
    guides(color=guide_legend(title="K-means cluster")) +
    theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) +
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

path = get_pwm_profile_plot_file_path(5)    
ggsave(paste0(path, '/profile.pdf'), lines, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)

lines2 = ggplot(cluster_data) +
    geom_line(aes(x = trim_length, y = pwm_score/log(10), group = gene), size = 3, alpha = 0.5) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Trim length') +
    ylab('-log10(PWM score)') +
    theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) +
    background_grid(major = 'xy') 

path = get_pwm_profile_plot_file_path(5)    
ggsave(paste0(path, '/profile_no_cluster.pdf'), lines2, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)

profile_abs_mag = cluster_data[, sum(abs(pwm_score/log(10))), by = gene]
setnames(profile_abs_mag, 'V1', 'pwm_score_abs_sum')
avg = mean(profile_abs_mag$pwm_score_abs_sum)

hist = ggplot(profile_abs_mag) +
    geom_histogram(aes(x = pwm_score_abs_sum), alpha = 0.7) +
    geom_vline(xintercept = avg, color = 'blue', size = 2) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Absolute sum of -log10(PWM score) by gene') +
    theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) +
    background_grid(major = 'xy') 

path = get_pwm_profile_plot_file_path(0)    
ggsave(paste0(path, '/profile_hist.pdf'), hist, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)


# # compare average trim lengths across groups
# mean_motif_data = motif_data[, mean(p_trim_given_gene), by = .(gene, trim_length)]
# setnames(mean_motif_data, 'V1', 'mean_p_trim_given_gene')
# wide_trim = mean_motif_data %>%
#     pivot_wider(names_from = 'trim_length', values_from = 'mean_p_trim_given_gene')%>%
#     as.data.table()

# cluster_plot_trim = determine_cluster_count_plot(wide_trim)
# cluster_data_trim = cluster_data(wide_trim, 5, 'mean_p_trim_given_gene')


