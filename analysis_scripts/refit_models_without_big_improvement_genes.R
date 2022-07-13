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

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[12])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

RESIDUAL_COMPARE_FEATURE <<- NULL
TYPE <<- 'full_v_gene_family_loss'
source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/residual_comparison_functions.R')
source('analysis_scripts/pwm_profile_functions.R')
source('scripts/model_evaluation_functions.R')

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


if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
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
if (grepl('_side_terminal', MODEL_TYPE1, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE1, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE1, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
    model1 = paste0(MODEL_TYPE1, '_', LEFT_SIDE1, '_length_melting_left')
} else {
    model1 = MODEL_TYPE1
}

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
    model2 = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
} else {
    model2 = MODEL_TYPE 
}

path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, paste0('compare_', model1, '-', model2) , paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
dir.create(path, recursive = TRUE)
 
# get difference
wide = together[, -c('p_gene')] %>% pivot_wider(names_from = 'model', values_from = 'rmse') %>% as.data.table()
slopes = wide[,(get(MODEL_TYPE) - get(MODEL_TYPE1)), by = gene]
setnames(slopes, 'V1', 'slope')
together = merge(together, slopes)

together = together[order(abs(slope))]

# get outliers
if (TRIM_TYPE == 'v_trim'){
    out_slope = -0.12
    outliers = together[slope < out_slope]
    no_outliers = together[slope > out_slope]
} else {
    outlier_slope_values = boxplot.stats(together[model %like% 'motif']$slope)$out
    outliers = together[slope %in% outlier_slope_values]
    no_outliers = together[!(slope %in% outlier_slope_values)]
}

# re-fit_model without outliers
motif_data = aggregate_all_subject_data()
new_motif_data = motif_data[gene %in% no_outliers$gene]

new_coefs = fit_model_by_group(new_motif_data, write_coeffs=FALSE)

# get coefficient plots 

dir.create(file.path(path, 'model_without_motif_addition_outliers'))
plot_melting_coefficient_heatmap_single_group(new_coefs, file_name = file.path(path, 'model_without_motif_addition_outliers', 'melting_heatmap.pdf'), limits = NULL, write_plot = TRUE, with_values = TRUE)
plot_distance_coefficient_heatmap_single_group(new_coefs, file_name = file.path(path, 'model_without_motif_addition_outliers', 'distance_heatmap.pdf'), limits = NULL, write_plot = TRUE, with_values = TRUE)
plot_model_coefficient_heatmap_single_group(new_coefs, file_name = file.path(path, 'model_without_motif_addition_outliers', 'motif_heatmap.pdf'), limits = NULL, write_plot = TRUE, with_values = TRUE)
plot_base_count_coefficient_heatmap_single_group(new_coefs, file_name = file.path(path, 'model_without_motif_addition_outliers', 'base_count_heatmap.pdf'), limits = NULL, write_plot = TRUE, with_values = TRUE)


# plot predicted distributions on held out (outlier genes)
model = fit_model(new_motif_data)
outlier_motif_data = motif_data[gene %in% outliers$gene]
outlier_motif_data = process_data_for_model_fit(outlier_motif_data)
outlier_motif_data$predicted_prob = temp_predict(model, newdata = outlier_motif_data)
outlier_motif_data[, empirical_prob := count/sum(count), by = .(subject, gene)]

dir.create(file.path(path, 'model_without_motif_addition_outliers', 'predicted_gene_dists'))
for (gene in unique(outlier_motif_data$gene)){
    plot_predicted_trimming_dists_single_group(outlier_motif_data, gene, file_name = file.path(path, 'model_without_motif_addition_outliers', 'predicted_gene_dists', paste0(gene, '.pdf')))
}

# per_gene_outlier_resid_motif_model = calculate_rmse_by_gene(outlier_motif_data) 
# per_gene_outlier_resid_NO_motif_model = calculate_rmse_by_gene(predicted_trims1[gene %in% outliers$gene]) 

# No fit without outliers and genes similar to outliers...
# re-fit_model without outliers
gene_families = get_gene_families(cluster_count = 4, combine_by_terminal = FALSE, full_sequence = TRUE, align = TRUE)$cluster_data
outlier_clusters = unique(gene_families[gene %in% unique(outliers$gene)]$clusters_grouped)
outlier_cluster_genes = unique(gene_families[clusters_grouped %in% outlier_clusters]$gene)
no_outliers = together[!(gene %in% outliers) & !(gene %in% outlier_cluster_genes)]
new_motif_data = motif_data[gene %in% no_outliers$gene]

new_coefs = fit_model_by_group(new_motif_data, write_coeffs=FALSE)

# get coefficient plots 

dir.create(file.path(path, 'model_without_motif_addition_outliers_and_families'))
plot_melting_coefficient_heatmap_single_group(new_coefs, file_name = file.path(path, 'model_without_motif_addition_outliers_and_families', 'melting_heatmap.pdf'), limits = NULL, write_plot = TRUE, with_values = TRUE)
plot_distance_coefficient_heatmap_single_group(new_coefs, file_name = file.path(path, 'model_without_motif_addition_outliers_and_families', 'distance_heatmap.pdf'), limits = NULL, write_plot = TRUE, with_values = TRUE)
plot_model_coefficient_heatmap_single_group(new_coefs, file_name = file.path(path, 'model_without_motif_addition_outliers_and_families', 'motif_heatmap.pdf'), limits = NULL, write_plot = TRUE, with_values = TRUE)
plot_base_count_coefficient_heatmap_single_group(new_coefs, file_name = file.path(path, 'model_without_motif_addition_outliers_and_families', 'base_count_heatmap.pdf'), limits = NULL, write_plot = TRUE, with_values = TRUE)


# plot predicted distributions on held out (outlier genes)
model = fit_model(new_motif_data)
outlier_motif_data = motif_data[gene %in% outliers$gene]
outlier_motif_data = process_data_for_model_fit(outlier_motif_data)
outlier_motif_data$predicted_prob = temp_predict(model, newdata = outlier_motif_data)
outlier_motif_data[, empirical_prob := count/sum(count), by = .(subject, gene)]

dir.create(file.path(path, 'model_without_motif_addition_outliers_and_families', 'predicted_gene_dists'))
for (gene in unique(outlier_motif_data$gene)){
    plot_predicted_trimming_dists_single_group(outlier_motif_data, gene, file_name = file.path(path, 'model_without_motif_addition_outliers_and_families', 'predicted_gene_dists', paste0(gene, '.pdf')))
}

