source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
library(cowplot)
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

ANNOTATION_TYPE <<- 'igor'
TRIM_TYPE <<- 'v_trim'
PRODUCTIVITY <<- 'nonproductive'
MOTIF_TYPE <<- 'unbounded' 
NCPU <<- 2
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
MODEL_GROUP <<- 'all_subjects'
GENE_WEIGHT_TYPE <<- 'p_gene_marginal' 
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 1
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 2
UPPER_TRIM_BOUND <<- 14
LOWER_TRIM_BOUND <<- 2
MODEL_TYPE <<- 'motif_two-side-base-count-beyond'
LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10

# load BCR data
VALIDATION_DATA_DIR <<- args[1]
VALIDATION_TYPE <<- 'validation_data_igh'
VALIDATION_TRIM_TYPE <<- args[2]
VALIDATION_PRODUCTIVITY <<- 'nonproductive'
VALIDATION_GENE_NAME <<- paste0(substring(VALIDATION_TRIM_TYPE, 1, 1), '_gene')
RESIDUAL_COMPARE_FEATURE <<- NULL

source('plotting_scripts/plotting_functions.R')
source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/residual_comparison_functions.R')
source('analysis_scripts/pwm_profile_functions.R')

if (MODEL_TYPE != 'null') {
    model = load_model()
} else {
    model = 'null'
} 
pwm = get_model_coefficient_data()

source(paste0('scripts/sampling_procedure_functions/p_gene_given_subject.R'), local = TRUE)
predicted_trims = get_predicted_distribution_data() 
predicted_trims = calculate_subject_gene_weight(predicted_trims)
training_pwm_scores = get_all_pwm_score(predicted_trims, pwm, NULL)
# training_avg = compare_weight_by_motif(predicted_trims[!(gene_count <= UPPER_TRIM_BOUND)])
training_avg = compare_weight_by_motif(predicted_trims)
training = merge(training_avg, unique(training_pwm_scores[, -c('trim_length')]), by = c('motif', 'gene'))

ANNOTATION_TYPE <<- VALIDATION_TYPE
TYPE <<- 'validation_data' 
TRIM_TYPE <<- VALIDATION_TRIM_TYPE
GENE_NAME <<- VALIDATION_GENE_NAME
PRODUCTIVITY <<- VALIDATION_PRODUCTIVITY

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('scripts/model_evaluation_functions.R')

validation_data = aggregate_validation_data(directory = VALIDATION_DATA_DIR)
validation_data$predicted_prob = temp_predict(model, newdata = validation_data)
validation_data[, empirical_prob := count/sum(count), by = .(subject, gene)]
source(paste0('scripts/sampling_procedure_functions/p_gene_given_subject.R'), local = TRUE)
validation_data = calculate_subject_gene_weight(validation_data)
validation_pwm_scores = get_all_pwm_score(validation_data, pwm, NULL)
# validation_avg = compare_weight_by_motif(validation_data[!(gene_count <= UPPER_TRIM_BOUND)])
validation_avg = compare_weight_by_motif(validation_data)
validation = merge(validation_avg, unique(validation_pwm_scores[, -c('trim_length')]), by = c('motif', 'gene'))

#get path
ANNOTATION_TYPE <<- 'igor'
path = get_manuscript_path()
path = file.path(path, 'motif_analysis')
dir.create(path)

together = rbind(training, validation)
together[gene %like% 'IGH', type := 'IGH']
together[gene %like% 'TRB', type := 'TRB']

# plot weight by motif

plot = ggplot(together) +
    geom_boxplot(aes(x = motif, y = total_weight_by_gene, fill = pwm_score))+
    geom_jitter(aes(x = motif, y = total_weight_by_gene), alpha = 0.3, width = 0.1)+
    facet_grid(rows = vars(type)) +
    scale_fill_distiller(palette = 'PuOr', name = 'PWM motif score', limits = c(-1.2, 1.2)) +
    xlab('Motif') +
    ylab('Total sequence weight (by gene)') +
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 25), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 16)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

ggsave(file.path(path, 'total_weight_by_motif.pdf'), plot = plot, width = 20, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

# plot total weight by PWM score
plot2 = ggplot(together) +
    geom_point(aes(x = total_weight_by_gene, y = pwm_score), alpha = 0.4, size = 3) +
    facet_grid(cols = vars(type)) +
    xlab('Total sequence weight (by gene)') +
    ylab('PWM motif score') +
    theme_cowplot(font_family = 'Arial') + 
    theme(legend.position = "none", text = element_text(size = 25), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), axis.text.x = element_text(size = 20)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

ggsave(file.path(path, 'total_weight_by_PWM_score.pdf'), plot = plot2, width = 14, height = 7, units = 'in', dpi = 750, device = cairo_pdf)

summed = together[, sum(total_weight_by_gene), by = .(motif, type, pwm_score)]
setnames(summed, 'V1', 'total_weight_all')
summed[type == 'IGH', total_weight_all := -1*total_weight_all]
ordered_weights = unique(summed[order(summed$pwm_score)]$motif)
ordered_weights = c(ordered_weights, '')
spacing = data.table(motif = c('', ''), type = c('TRB', 'IGH'), pwm_score = c(NA, NA), total_weight_all = c(0.35, -0.35))
summed = rbind(summed, spacing)
summed$motif = factor(summed$motif, levels = ordered_weights)
require(ggpol)
plot3 = ggplot(summed, aes(y = total_weight_all, x = motif, fill = pwm_score)) +
    geom_bar(position = position_dodge(width=1), stat='identity') + 
    facet_share(~type, dir = "h", scales = "free", reverse_num = TRUE) +
    coord_flip() +
    xlab('') +
    ylab('Motif frequency') +
    scale_fill_distiller(palette = 'PuOr', name = 'PWM motif score', limits = c(-1.25, 1.25), na.value="transparent") +
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 25), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), axis.text.x = element_text(size = 20)) + 
    guides(fill = guide_colourbar(barwidth = 2, barheight = 15))+
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

ggsave(file.path(path, 'total_weight_barplot.pdf'), plot = plot3, width = 16, height = 20, units = 'in', dpi = 750, device = cairo_pdf)

summed2 = summed %>% pivot_wider(names_from = 'type', values_from = 'total_weight_all') %>% as.data.table()
summed2$IGH = -1*(summed2$IGH)
summed2[is.na(TRB), TRB := 0]
summed2[is.na(IGH), IGH := 0]

plot4 = ggplot(summed2, aes(x = TRB, y = IGH, color = pwm_score)) +
    geom_point(alpha = 0.8, size = 5) +
    xlab('TRB motif frequency') +
    ylab('IGH motif frequency') +
    geom_abline(intercept = 0, size = 3, color = 'black') +
    scale_color_distiller(palette = 'PuOr', name = 'PWM motif score', limits = c(-1.25, 1.25), na.value="transparent") +
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 25), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), axis.text.x = element_text(size = 20)) + 
    guides(color = guide_colourbar(barwidth = 2, barheight = 15))+
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

ggsave(file.path(path, 'total_weight_scatter.pdf'), plot = plot4, width = 14, height = 12, units = 'in', dpi = 750, device = cairo_pdf)

valid_count = compare_weight_by_base_count(validation_data)
valid_count_pwm = get_all_base_count_pwm_score(validation_data, pwm)
valid = merge(valid_count, valid_count_pwm)
valid$type = 'IGH'
train_count = compare_weight_by_base_count(predicted_trims)
train_count_pwm = get_all_base_count_pwm_score(predicted_trims, pwm)
train = merge(train_count, train_count_pwm)
train$type = 'TRB'
together2 = rbind(valid, train)
together2[type == 'IGH', total_weight := -1* total_weight]
together2[, base_counts := paste0(left_base_count_GC, ', ', right_base_count_GC, ', ', right_base_count_AT)]
ordered_weights = unique(together2[order(together2$pwm_score)]$base_counts)
ordered_weights = c(ordered_weights, '')
spacing = data.table(base_counts = c('', ''), type = c('TRB', 'IGH'), pwm_score = c(NA, NA), total_weight = c(0.39, -0.39))
together2 = rbind(together2, spacing, fill = TRUE)
together2$base_counts = factor(together2$base_counts, levels = ordered_weights)

plot5 = ggplot(together2, aes(y = total_weight, x = base_counts, fill = pwm_score)) +
    geom_bar(position = position_dodge(width=1), stat='identity') + 
    facet_share(~type, dir = "h", scales = "free", reverse_num = TRUE) +
    coord_flip() +
    xlab('Base counts: 5\'GC, 3\'GC, 3\'AT') +
    ylab('Base count frequency') +
    scale_fill_distiller(palette = 'PuOr', name = 'PWM motif score', limits = c(-2.64, 2.64), na.value="transparent") +
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 25), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), axis.text.x = element_text(size = 20)) + 
    guides(fill = guide_colourbar(barwidth = 2, barheight = 15))+
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

ggsave(file.path(path, 'total_weight_barplot_base_counts.pdf'), plot = plot5, width = 16, height = 45, units = 'in', dpi = 750, device = cairo_pdf)


