source('config/config.R')

library(grid)
library(cowplot)
library(ggplot2)
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

# NOTE: This method is only applicable for models fit across all subjects!
MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- args[6]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[7])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[8])

UPPER_TRIM_BOUND <<- as.numeric(args[9]) 

MODEL_TYPE <<- 'motif_distance_two_side_terminal_melting' 

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[10])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('analysis_scripts/bootstrap_analysis_functions.R')
source('plotting_scripts/plotting_functions.R')

motif_data = aggregate_all_subject_data()

together = data.table()
for (bootstrap in seq(1, 99)){
    pwm_file_path = get_pwm_matrix_file_path_bootstrap(bootstrap)
    pwm_file_name = get_pwm_matrix_file_name(subgroup = NULL)
    location = file.path(pwm_file_path, pwm_file_name)
    data = fread(location) 
    long = data %>% pivot_longer(cols = starts_with("motif"), names_to = 'position', values_to = 'coef') %>% as.data.table()
    long$bootstrap_dataset = bootstrap
    together = rbind(together, long)
}

original = get_model_coefficient_data() 
original = original[parameter %like% 'motif']
setnames(original, 'parameter', 'position')
setnames(original, 'coefficient', 'coef')
together = rbind(together, original, fill = TRUE)

position_values = map_positions_to_values(unique(together$position))
together = merge(together, position_values, by.x = 'position', by.y = 'positions')

# convert to log_10
together$log_10_pdel = together$coef/log(10)
# order variables
together$base = factor(together$base, levels = c('T', 'G', 'C', 'A'))
together$values = factor(together$values)
together = data.table(together)

plot = ggplot(together[!is.na(bootstrap_dataset)], aes(x=log_10_pdel)) +
    facet_grid(rows = vars(base), cols = vars(values)) +
    geom_histogram(aes(y = stat(width*density), group = values), position = 'identity')+
    geom_vline(data = together[is.na(bootstrap_dataset)], aes(xintercept=log_10_pdel), color = '#1b9e77', linetype = 'F1', size=1) +
    theme_cowplot(font_family = 'Arial') + 
    xlab('log10(probability of deletion)') +
    ylab ('Proportion of bootstrapped datasets') +
    theme(text = element_text(size = 20)) +
    background_grid(major = 'xy') +
    panel_border(color = 'gray60', size = 1.5)

path = get_coef_variations_file_path()
file = 'coef_variations_bootstrapped_datasets.pdf'

file_name = file.path(path, file)
ggsave(file_name, plot = plot, width = 15, height = 10, units = 'in', dpi = 750, device = cairo_pdf)


motif_length = RIGHT_NUC_MOTIF_COUNT + LEFT_NUC_MOTIF_COUNT

groups = list(seq(1, 20), seq(21, 40), seq(41, 60), seq(61, 80), seq(81, 99)) 
for (group in groups){
    temp_data = together[!is.na(bootstrap_dataset) & bootstrap_dataset %in% group]
    temp_data$base = factor(temp_data$base, levels = c('T', 'G', 'C', 'A'))
    plot = ggplot(temp_data, aes(x=values, y=base, fill=log_10_pdel)) +
        facet_wrap(~bootstrap_dataset, ncol = 5, nrow = 4) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Position') +
        ylab ('Base') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank(), panel.spacing.y = unit(1, "lines"), panel.spacing.x = unit(2, "lines")) +
        geom_vline(xintercept = LEFT_NUC_MOTIF_COUNT + 0.5, size = 3, color = 'black') +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)') +
        coord_cartesian(ylim = c(1, 4), clip = "off")

    path = get_coef_heatmap_file_path()
    file_name = paste0(path, '/bootstrapped_genes_heatmaps_', group[1], '.pdf')

    ggsave(file_name, plot = plot, width = 18, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
}

vardata = together[!is.na(bootstrap_dataset), var(log_10_pdel), by = .(position, base, values)]
vardata$base = factor(vardata$base, levels = c('T', 'G', 'C', 'A'))
setnames(vardata, 'V1', 'var_log_10_pdel')
plot = ggplot(vardata, aes(x=values, y=base, fill=var_log_10_pdel)) +
    geom_tile() +
    theme_cowplot(font_family = 'Arial') + 
    xlab('Position') +
    ylab ('Base') +
    theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
    geom_vline(xintercept = LEFT_NUC_MOTIF_COUNT + 0.5, size = 3, color = 'black') +
    guides(fill = guide_colourbar(barheight = 14)) +
    scale_fill_viridis_c(name = 'Variance of\nlog10(probability of deletion)') +
    coord_cartesian(ylim = c(1, 4), clip = "off")

path = get_coef_heatmap_file_path()
file_name = paste0(path, '/bootstrapped_genes_variance_heatmap.pdf')

ggsave(file_name, plot = plot, width = 9, height = 4.1, units = 'in', dpi = 750, device = cairo_pdf)


# plot small multiples heatmaps of background
backgrounds = data.table()
for (bootstrap in seq(1, 99)){
    gene_path = get_data_genes_path_bootstrap(bootstrap)
    genes = fread(paste0(gene_path, '/sampled_genes.tsv')) 
    boot_data = get_bootstrap_sample_from_gene_list(motif_data, genes) 
    boot_background = get_base_composition_counts(boot_data)
    long = boot_background %>% pivot_longer(cols = starts_with("motif"), names_to = 'position', values_to = 'freq') %>% as.data.table()
    long$bootstrap_dataset = bootstrap
    backgrounds = rbind(backgrounds, long)
}

position_values = map_positions_to_values(unique(backgrounds$position))
backgrounds_together = merge(backgrounds, position_values, by.x = 'position', by.y = 'positions')

# order variables
backgrounds_together$base = factor(backgrounds_together$base, levels = c('T', 'G', 'C', 'A'))
backgrounds_together$values = factor(backgrounds_together$values)

# get log fold change above 0.25 background freq
backgrounds_together$log_fold_freq = log2(backgrounds_together$freq) - log2(0.25)
motif_length = RIGHT_NUC_MOTIF_COUNT + LEFT_NUC_MOTIF_COUNT

max = round(max(abs(backgrounds_together$log_fold_freq)), 1) + 0.1

groups = list(seq(1, 20), seq(21, 40), seq(41, 60), seq(61, 80), seq(81, 99)) 
for (group in groups){
    plot = ggplot(backgrounds_together[bootstrap_dataset %in% group], aes(x=values, y=base, fill=log_fold_freq)) +
        facet_wrap(~bootstrap_dataset, ncol = 5, nrow = 4)+
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Position') +
        ylab ('Base') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        geom_vline(xintercept = LEFT_NUC_MOTIF_COUNT + 0.5, size = 3, color = 'black') +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log2 fold change in\nbackground base frequency', limits = c(-1*max, max)) +
        coord_cartesian(ylim = c(1, 4), clip = "off")
    
    path = get_coef_heatmap_file_path()
    file_name = paste0(path, '/background_bases_bootstrapped_genes_heatmaps_', group[1], '.pdf')

    ggsave(file_name, plot = plot, width = 18, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
} 
