source('config/config.R')

library(cli)
library(devtools)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(foreach)
library(doParallel)
library(tidyverse)
library(plyr)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)


ANNOTATION_TYPE <<- 'igor' 
TRIM_TYPE <<- 'v_trim'
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
PRODUCTIVITY <<- 'nonproductive'

MOTIF_TYPE <<- 'unbounded' 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(2)

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- 'p_gene_marginal' 
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

UPPER_TRIM_BOUND <<- as.numeric(14) 
LOWER_TRIM_BOUND <<- 2 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(10)

TYPE <<- 'log_loss'
LOSS_GENE_WEIGHT <<- 'p_gene_given_subject' 
stopifnot(LOSS_GENE_WEIGHT %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform', 'p_gene_marginal_all_seqs'))

source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/model_evaluation_functions.R')
source('analysis_scripts/rel_importance_functions.R')

all_eval_results = data.table()
rel_import_results = data.table()
# annotation_types = c('validation_data_alpha', 'validation_data_beta', 'validation_data_gamma')
annotation_types = c('validation_data_alpha', 'validation_data_beta', 'igor', 'validation_data_gamma', 'validation_data_delta', 'validation_data_igh')

trim_types = c('j_trim', 'v_trim')
prods = c('productive', 'nonproductive')

for (trim_type in trim_types){
    for (type in annotation_types) {
        for (prod in prods){
            if (prod == 'nonproductive' & type %in% c('validation_data_gamma', 'validation_data_delta')) {
                next
            }
            if (prod == 'productive' & type == 'validation_data_igh'){
                next
            }
            VALIDATION_TYPE <<- type
            VALIDATION_PRODUCTIVITY <<- prod
            VALIDATION_TRIM_TYPE <<- trim_type

            source('analysis_scripts/rel_importance_functions.R')

            temp_rel_results = compile_rel_importance_results()
            temp_rel_results$loss_type = type
            temp_rel_results$trim_type = trim_type
            temp_rel_results$productivity = prod
            rel_import_results = rbind(rel_import_results, temp_rel_results, fill = TRUE)
        }
    }
}


for (trim_type in trim_types){
    for (type in annotation_types) {
        for (prod in prods){
            if (prod == 'nonproductive' & type %in% c('validation_data_gamma', 'validation_data_delta')) {
                next
            }
            if (prod == 'productive' & type == 'validation_data_igh'){
                next
            }
            ANNOTATION_TYPE <<- type
            TRIM_TYPE <<- trim_type
            GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
            PRODUCTIVITY <<- prod

            source('scripts/model_evaluation_functions.R')
            source('plotting_scripts/plotting_functions.R')
            source('plotting_scripts/model_evaluation_functions.R')
            temp_eval_results = compile_evaluation_results(type)
            setnames(temp_eval_results, type, 'loss')
            temp_eval_results$loss_type = type
            temp_eval_results$trim_type = trim_type
            temp_eval_results$productivity = prod
            all_eval_results = rbind(all_eval_results, temp_eval_results, fill = TRUE)
        }
    }
}

together = merge(all_eval_results[model_type == 'motif_two-side-base-count-beyond'], rel_import_results, by = c('loss_type', 'trim_type', 'productivity'))
# load relative importance data

plot = ggplot(together) +
    geom_abline(intercept = 0, size = 3, color = 'black')+ 
    geom_point(aes(x = base_count_score, y = motif_score, color = loss), size = 7) +
    facet_grid(rows = vars(trim_type), cols = vars(productivity))+
    theme_cowplot(font_family = 'Arial') + 
    xlab('\nBase-count-beyond feature scale') +
    ylab('Motif feature scale\n')+ 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.1) +
    theme(text = element_text(size = 38), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 38), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), strip.text = element_text(size = 38), panel.spacing = unit(3, "lines"), legend.key.height = unit(4, 'cm'), legend.key.width = unit(1.5, 'cm'))  +
    scale_color_viridis_c(name = 'Expected\nper-sequence\nlog loss') 

ANNOTATION_TYPE <<- 'igor'
path = get_manuscript_path()
file_name = paste0(path, '/feature_scale_comparison.pdf')
ggsave(file_name, plot = plot, width = 25, height = 20, units = 'in', dpi = 750, device = cairo_pdf)

together[, ratio := base_count_score/motif_score]
plot2 = ggplot(together) + 
    geom_point(aes(x = ratio, y = loss, shape = trim_type, color = loss_type), size = 6)+
    facet_grid(cols = vars(productivity))+
    theme_cowplot(font_family = 'Arial') + 
    xlab('\nFeature scale ratio (base count/motif)') +
    ylab('Expected per-sequence log loss\n')+ 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.1) +
    theme(text = element_text(size = 38), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 38), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), strip.text = element_text(size = 38), panel.spacing = unit(3, "lines"))

file_name2 = paste0(path, '/feature_scale_scatter.pdf')
ggsave(file_name2, plot = plot2, width = 20, height = 9, units = 'in', dpi = 750, device = cairo_pdf)

plot3 = ggplot(together) +
    geom_abline(intercept = 0, size = 3, color = 'black')+ 
    geom_point(aes(x = base_count_score, y = motif_score, color = loss_type), size = 7) +
    facet_grid(rows = vars(trim_type), cols = vars(productivity))+
    theme_cowplot(font_family = 'Arial') + 
    xlab('\nBase-count-beyond feature scale') +
    ylab('Motif feature scale\n')+ 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.1) +
    theme(text = element_text(size = 38), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 38), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), strip.text = element_text(size = 38), panel.spacing = unit(3, "lines"), legend.key.height = unit(4, 'cm'), legend.key.width = unit(1.5, 'cm')) 

file_name3 = paste0(path, '/feature_scale_comparison_types.pdf')
ggsave(file_name3, plot = plot3, width = 25, height = 20, units = 'in', dpi = 750, device = cairo_pdf)
