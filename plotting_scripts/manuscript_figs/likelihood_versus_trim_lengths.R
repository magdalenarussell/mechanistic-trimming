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

MODEL_TYPE <<- 'motif_two-side-base-count-beyond'

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

UPPER_TRIM_BOUND <<- as.numeric(14) 
LOWER_TRIM_BOUND <<- 2 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(10)

TYPE <<- 'log_loss'
RESIDUAL_COMPARE_FEATURE <<- NULL

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/model_evaluation_functions.R')
source('plotting_scripts/residual_comparison_functions.R')

if (MODEL_TYPE != 'null') {
    model = load_model()
} else {
    model = 'null'
} 

# Compile data for all subjects
igor_nonproductive_data = aggregate_all_subject_data()

source('config/validation_file_paths.R')

seq_weights = data.table()
log_lik = data.table()
model_prob = data.table()
gene_use = data.table()
rmse_all = data.table()
for (valid in c('igor', 'validation_data_igh', 'validation_data_delta', 'validation_data_gamma')){
    for (type in c('productive', 'nonproductive')){
        if (valid == 'igor' & type == 'nonproductive'){
            assign(paste0(valid, '_', type, '_data'), igor_nonproductive_data)
        } else if (valid == 'validation_data_delta' & type == 'nonproductive'){
            next
        } else if (valid == 'validation_data_gamma' & type == 'nonproductive'){
            next
        } else {
            IGOR = TCR_REPERTOIRE_DATA_igor
            if (valid == 'validation_data_igh'){
                VALIDATION_DATA_DIR <<- get(paste0(toupper(valid), '_', type))
            } else {
                VALIDATION_DATA_DIR <<- get(toupper(valid))
            }
            VALIDATION_TYPE <<- valid
            VALIDATION_TRIM_TYPE <<- 'v_trim'
            VALIDATION_PRODUCTIVITY <<- type
            VALIDATION_GENE_NAME <<- paste0(substring(VALIDATION_TRIM_TYPE, 1, 1), '_gene')

            ANNOTATION_TYPE <<- VALIDATION_TYPE
            TYPE <<- 'validation_data' 
            TRIM_TYPE <<- VALIDATION_TRIM_TYPE
            GENE_NAME <<- VALIDATION_GENE_NAME
            PRODUCTIVITY <<- VALIDATION_PRODUCTIVITY

            source('scripts/data_compilation_functions.R')
            source('scripts/model_fitting_functions.R')
            source('scripts/model_evaluation_functions.R')

            assign(paste0(valid, '_', type, '_data'), aggregate_validation_data(directory = VALIDATION_DATA_DIR))
        } 
        source(paste0('scripts/sampling_procedure_functions/p_gene_marginal.R'), local = TRUE)
        data0 = get(paste0(valid, '_', type, '_data'))
        data0 = calculate_subject_gene_weight(data0)

        source(paste0('scripts/sampling_procedure_functions/p_gene_given_subject.R'), local = TRUE)
        data = get(paste0(valid, '_', type, '_data'))
        data = calculate_subject_gene_weight(data)

        if (MODEL_TYPE != 'null') {
            data$prediction = temp_predict(model, newdata = data)
            data0$prediction = temp_predict(model, newdata = data0)
        } else {
            data$prediction = exp(1)/sum(rep(exp(1), UPPER_TRIM_BOUND -LOWER_TRIM_BOUND + 1))
            data0$prediction = exp(1)/sum(rep(exp(1), UPPER_TRIM_BOUND -LOWER_TRIM_BOUND + 1))
        }

        setnames(data0, 'prediction', 'predicted_prob')
        setnames(data0, 'p_trim_given_gene', 'empirical_prob')
        rmse = calculate_rmse_by_gene(data0) 
        rmse$type = type
        rmse$validation_type = valid
        rmse_all = rbind(rmse_all, rmse)

        data = data[count != 0]
        data[, log_pred := log(prediction)]
        data[, log_lik := log_pred * count]
        data[, weighted_log_lik := log_pred * weighted_observation]
        data[, total_gene_count := sum(count), by = .(gene)]
        data[, type := type]
        data[, validation_type := valid] 
        # average across people
        total_subj = length(unique(data$subject))
        temp = data[, sum(weighted_observation), by = .(trim_length, gene, type, validation_type)]
        setnames(temp, 'V1', 'total_weighted_obs')
        seq_weights = rbind(seq_weights, temp)
        temp2 = data[, sum(weighted_log_lik), by = .(trim_length, gene, type, validation_type)]
        setnames(temp2, 'V1', 'weighted_log_lik')
        log_lik = rbind(log_lik, temp2)
        cols = c('log_pred', 'gene', 'trim_length', 'type', 'validation_type')
        temp3 = unique(data[, ..cols])
        model_prob = rbind(model_prob, temp3)
        temp4 = data[, sum(p_gene), by = .(gene, type, validation_type)]
        temp4$p_gene_marg = temp4$V1/total_subj
        gene_use = rbind(gene_use, temp4)
    }
}

cols = c('gene', 'rmse', 'type', 'validation_type')
rmse_all = unique(rmse_all[, ..cols])
gene_use = merge(gene_use, rmse_all[type == 'productive'], by = c('gene', 'validation_type'))
gene_use_np = gene_use[type.x == 'nonproductive']
gene_use_p = gene_use[type.x == 'productive']
gene_use = merge(gene_use_np, gene_use_p)

#get path
PRODUCTIVITY <<- 'nonproductive'
ANNOTATION_TYPE <<- 'igor'
path = get_manuscript_path()
path = file.path(path, 'prod_nonprod_analysis')
dir.create(path)

plot = ggplot(seq_weights, aes(x = trim_length, y =total_weighted_obs)) +
    geom_point(alpha = 0.5, size = 3) +
    facet_grid(rows = vars(type), cols = vars(validation_type))+
    # scale_y_continuous(trans = scales::log_trans(), 
                       # breaks = scales::log_breaks())+
    geom_smooth(method = 'lm', size = 2) +
    xlab('Trim length') +
    ylab('Total sequence weight (by gene/length)') +
    theme_cowplot(font_family = 'Arial') + 
    theme(legend.position = "none", text = element_text(size = 25), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), axis.text.x = element_text(size = 20)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

ggsave(file.path(path, 'avg_weight_by_trim_length.pdf'), plot = plot, width = 40, height = 15, units = 'in', dpi = 750, device = cairo_pdf)

plot2 = ggplot(log_lik, aes(x = trim_length, y =weighted_log_lik)) +
    geom_point(alpha = 0.5, size = 3) +
    facet_grid(rows = vars(type), cols = vars(validation_type))+
    # scale_y_continuous(trans = scales::log_trans(), 
                       # breaks = scales::log_breaks())+
    geom_smooth(method = 'lm', size = 2) +
    xlab('Trim length') +
    ylab('Total weighted log likelihood (by gene/length)') +
    theme_cowplot(font_family = 'Arial') + 
    theme(legend.position = "none", text = element_text(size = 25), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), axis.text.x = element_text(size = 20)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

ggsave(file.path(path, 'avg_weighted_log_lik_by_trim_length.pdf'), plot = plot2, width = 40, height = 15, units = 'in', dpi = 750, device = cairo_pdf)

plot3 = ggplot(model_prob, aes(x = trim_length, y = log_pred)) +
    geom_point(alpha = 0.5, size = 3) +
    facet_grid(rows = vars(type), cols = vars(validation_type))+
    # scale_y_continuous(trans = scales::log_trans(), 
                       # breaks = scales::log_breaks())+
    xlab('Trim length') +
    ylab('Model log probability') +
    theme_cowplot(font_family = 'Arial') + 
    theme(legend.position = "none", text = element_text(size = 25), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), axis.text.x = element_text(size = 20)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

ggsave(file.path(path, 'pred_by_trim_length.pdf'), plot = plot3, width = 40, height = 15, units = 'in', dpi = 750, device = cairo_pdf)

plot4 = ggplot(gene_use, aes(x = log(p_gene_marg.x), y = log(p_gene_marg.y), color = rmse.x)) +
    facet_grid(cols = vars(validation_type))+
    geom_point(alpha = 0.5, size = 4) +
    geom_abline(intercept = 0, size = 2) +
    xlab('Nonproductive gene usage (log)') +
    ylab('Productive gene usage (log)') +
    scale_color_viridis_c() +
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 25), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5), axis.text.x = element_text(size = 20)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

ggsave(file.path(path, 'gene_usage_rmse.pdf'), plot = plot4, width = 45, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

