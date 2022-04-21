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

# predict trimming probabilities using only pwm
predicted = predict_trimmming_given_pwm_scores(condensed)

together = merge(motif_data, predicted, by = c('gene', 'trim_length', 'motif', positions))
setnames(together, 'p_trim_given_gene', 'empirical_prob')
setnames(together, 'predicted_p_trim_given_gene', 'predicted_prob')

path = get_pwm_prediction_residual_plot_file_path()

file_name = paste0(path, '/all_residuals_outliers_colored.pdf')

# determine outliers
together$residual = together$empirical_prob - together$predicted_prob
mean_data = together[, mean(abs(residual)), by = .(gene)]
setnames(mean_data, 'V1', 'mean_residual')
outliers = boxplot.stats(mean_data$mean_residual)$out
out_ind = which(mean_data$mean_residual %in% c(outliers))
outlier_genes = mean_data[out_ind,]$gene

# get residual by trim and add outlier column
mean_data_trim = together[, mean(residual), by = .(gene, trim_length)]
setnames(mean_data_trim, 'V1', 'mean_residual')
mean_data_trim[, outlier := FALSE]
mean_data_trim[gene %in% outlier_genes, outlier := TRUE]

plot = ggplot() +
        geom_line(data = mean_data_trim, aes(x = trim_length, y = mean_residual, group = gene, color = outlier), size = 2, alpha = 0.7) +
        geom_hline(yintercept = 0, color = 'gray', size = 3) +
        ylab('Observed prob - Predicted prob\n') +
        theme_cowplot(font_family = 'Arial') + 
        theme(text = element_text(size = 30), axis.text.x=element_blank(), axis.title.x = element_blank(), axis.text.y = element_text(size = 20), axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        ylim(-0.5, 1)
 
ggsave(file_name, plot, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)

# add outlier column to condensed data
condensed[, outlier := FALSE]
condensed[gene %in% outlier_genes, outlier := TRUE]

profiles = ggplot(condensed) +
    geom_line(aes(x = trim_length, y = pwm_score/log(10), group = gene, color = outlier), size = 3, alpha = 0.5) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Trim length') +
    ylab('-log10(PWM score)') +
    theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) +
    background_grid(major = 'xy') 

path = get_pwm_profile_plot_file_path(0)    
ggsave(paste0(path, '/profile_no_cluster_outlier_colored.pdf'), profiles, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)


