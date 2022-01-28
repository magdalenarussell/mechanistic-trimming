source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(ggplot2)
library(cowplot)
library(ggpubr)
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

if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[12])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

RESIDUAL_COMPARE_FEATURE <<- NULL

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/residual_comparison_functions.R')
source('analysis_scripts/pwm_profile_functions.R')

# Read in dist data
predicted_trims1 = get_predicted_distribution_data() 
per_gene_resid1 = calculate_rmse_by_gene(predicted_trims1) 
per_gene_resid1$model = MODEL_TYPE

MODEL_TYPE1 = MODEL_TYPE
LEFT_SIDE1 = LEFT_SIDE_TERMINAL_MELT_LENGTH
MODEL_TYPE <<- args[12]
stopifnot(MODEL_TYPE %like% 'motif')
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[13])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[14])


if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[15])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')
source('analysis_scripts/pwm_profile_functions.R')

# Read in dist data
predicted_trims2 = get_predicted_distribution_data() 
per_gene_resid2 = calculate_rmse_by_gene(predicted_trims2) 
per_gene_resid2$model = MODEL_TYPE

#merge two predictions
together = rbind(per_gene_resid1, per_gene_resid2)

#get plot path and name
if (grepl('_side_terminal_melting', MODEL_TYPE1, fixed = TRUE)){
    model1 = paste0(MODEL_TYPE1, '_', LEFT_SIDE1, '_length_melting_left')
} else {
    model1 = MODEL_TYPE1
}

if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
    model2 = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
} else {
    model2 = MODEL_TYPE 
}

path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, paste0('compare_', model1, '-', model2) , paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
dir.create(path, recursive = TRUE)
    
plot = ggplot(together) +
    geom_point(aes(x = model, y = rmse, color = subject_count), size =6, alpha = 0.6) +
    geom_line(aes(x = model, y = rmse, group = gene, color = subject_count), size = 4, alpha = 0.6) +
    theme_cowplot(font_family = 'Arial') + 
    theme(legend.key.height = unit(3, 'cm'), text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

file_name = file.path(path, 'compare_rmse_by_subject_count.pdf') 
ggsave(file_name, plot = plot, width = 18, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

wide = together %>% pivot_wider(names_from = 'model', values_from = 'rmse') %>% as.data.table()
slopes = wide[,(get(MODEL_TYPE) - get(MODEL_TYPE1)), by = gene]
setnames(slopes, 'V1', 'slope')
together = merge(together, slopes)

together = together[order(abs(slope))]

plot2 = ggplot(together) +
    geom_line(aes(x = model, y = rmse, group = gene, color = slope), size = 4, alpha = 0.8) +
    geom_point(aes(x = model, y = rmse, color = slope), size =6, alpha = 0.8) +
    scale_color_gradient2() +
    theme_cowplot(font_family = 'Arial') + 
    theme(legend.key.height = unit(3, 'cm'), text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

file_name = file.path(path, 'compare_rmse_by_slope.pdf') 
ggsave(file_name, plot = plot2, width = 18, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

# Compile data for all subjects
positions = get_positions()
cols = c('gene', 'trim_length', 'motif', positions)
condensed = unique(predicted_trims2[,..cols]) 

# Read in model coefficient data 
pwm = get_model_coefficient_data() 

# calculate pwm score by gene and trim length
condensed[, pwm_score := unlist(lapply(motif, function(x) as.numeric(get_pwm_score(pwm, x, positions))))]

total_score = condensed[, mean(abs(pwm_score)), by = gene]
setnames(total_score, 'V1', 'mean_abs_pwm_score')

tog_score = merge(slopes, total_score)

plot3 = ggplot(tog_score) +
    geom_point(aes(x = slope, y = mean_abs_pwm_score), size = 6, alpha = 0.6)+
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
    background_grid(major = 'xy') + 
    ylab('Mean absolute value of PWM score') +
    panel_border(color = 'gray60', size = 1.5) 

file_name = file.path(path, 'compare_pwm_only_score_slope.pdf') 
ggsave(file_name, plot = plot3, width = 18, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

# NOW, compare slope to the pwm-only trimming prediction rmse
predicted = predict_trimmming_given_pwm_scores(condensed)
setnames(predicted, 'predicted_p_trim_given_gene', 'predicted_prob')

pwm_only_predictions = merge(predicted_trims2[, -c('predicted_prob')], predicted)
pwm_only_rmse = calculate_rmse_by_gene(pwm_only_predictions)

tog_rmse = merge(slopes, pwm_only_rmse)

plot4 = ggplot(tog_rmse) +
    geom_point(aes(x = slope, y = rmse), size = 6, alpha = 0.6)+
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
    background_grid(major = 'xy') + 
    ylab('PWM-only prediction RMSE') +
    panel_border(color = 'gray60', size = 1.5) 

file_name = file.path(path, 'compare_pwm_only_rmse_slope.pdf') 
ggsave(file_name, plot = plot4, width = 18, height = 10, units = 'in', dpi = 750, device = cairo_pdf)


