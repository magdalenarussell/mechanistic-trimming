source('mechanistic-trimming/config/config.R')

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

source(paste0(MOD_PROJECT_PATH,'plotting_scripts/plotting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'plotting_scripts/model_evaluation_functions.R'))

all_eval_results = data.table()

# assign all validation data, trim-type, and productivity variables
annotation_types = c('igor', 'validation_data_alpha', 'validation_data_beta', 'validation_data_gamma', 'validation_data_igh')
trim_types = c('j_trim', 'v_trim')
prods = c('productive', 'nonproductive')

# get all validation results
for (trim_type in trim_types){
    for (type in annotation_types) {
        for (prod in prods){
            ANNOTATION_TYPE <<- type
            TRIM_TYPE <<- trim_type
            GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
            PRODUCTIVITY <<- prod
            
            source('scripts/model_evaluation_functions.R')
            source('plotting_scripts/plotting_functions.R')
            source('plotting_scripts/model_evaluation_functions.R')

            temp_eval_results = compile_evaluation_results(type)
            setnames(temp_eval_results, type, 'loss', skip_absent = TRUE)
            temp_eval_results$loss_type = type
            temp_eval_results$trim_type = trim_type
            temp_eval_results$productivity = prod
            all_eval_results = rbind(all_eval_results, temp_eval_results, fill = TRUE)
        }
    }
}

# get model types, and make them neat
orig_model_types = c("motif_two-side-base-count-beyond", "motif", 'null')
new_model_types = c('1x2motif + two-side\nbase-count beyond\n(12 params)', '2x4motif (18 params)', 'null (0 params)')
all_eval_results = all_eval_results[model_type %in% orig_model_types]
all_eval_results$model_type = mapvalues(all_eval_results$model_type, from = orig_model_types, to = new_model_types) 
all_eval_results[loss_type == 'igor', loss_type := paste0(loss_type, '_', trim_type, '_', productivity)]
neat_names = make_model_names_neat(new_model_types)

# set plot colors
colors = set_color_palette(c(neat_names, '2x4motif (18 params)'), with_params = TRUE)

# make all loss type names neat and ordered
loss_types = c('igor_j_trim_productive', 'igor_v_trim_nonproductive', 'validation_data_alpha', 'validation_data_beta', 'validation_data_gamma', 'validation_data_igh', 'igor_j_trim_nonproductive', 'igor_v_trim_productive') 

nice_loss_types = c('TCRB training\ndataset J-genes\n(not included during\nmodel training)', 'full TCRB V-gene\ntraining dataset', 'TCRA testing\ndataset', 'TCRB testing\ndataset', 'TCRG testing\ndataset', 'IGH testing\ndataset', 'TCRB training\ndataset J-genes\n(not included during\nmodel training)', 'TCRB training\ndataset V-genes\n(not included during\nmodel training)')

loss_order = c('TCRB training\ndataset J-genes\n(not included during\nmodel training)', 'full TCRB V-gene\ntraining dataset', 'TCRB training\ndataset V-genes\n(not included during\nmodel training)', 'TCRB testing\ndataset', 'TCRA testing\ndataset', 'TCRG testing\ndataset', 'IGH testing\ndataset')

all_eval_results$nice_loss_type = mapvalues(all_eval_results$loss_type, from = loss_types, to=nice_loss_types)
all_eval_results[nice_loss_type == 'full TCRB V-gene\ntraining dataset', trim_type := 'v_trim']

# make trimming types neat and ordered
all_eval_results$trim_type = mapvalues(all_eval_results$trim_type, from = unique(all_eval_results$trim_type), to=c('J-gene trimming', 'V-gene trimming'))
all_eval_results$trim_type = factor(all_eval_results$trim_type, levels = c('V-gene trimming', 'J-gene trimming'))

# reformat model names to be nicer
all_eval_results = all_eval_results[motif_type == MOTIF_TYPE]
nice_names = make_model_names_neat(unique(all_eval_results$model_type)) 
all_eval_results$nice_model_type = mapvalues(all_eval_results$model_type, from = unique(all_eval_results$model_type), to = nice_names)

# order losses
ordered_losses = loss_order
all_eval_results$nice_loss_type = factor(all_eval_results$nice_loss_type, levels = ordered_losses)
last_loss = ordered_losses[length(ordered_losses)]
label_data = all_eval_results[nice_loss_type == last_loss] 

# save data for training data only (for reference)
motif_base_count_training_loss = all_eval_results[nice_loss_type == 'full TCRB V-gene\ntraining dataset' & nice_model_type == '1x2motif + two-side\nbase-count beyond\n(12 params)']$loss
motif_base_count_color = colors[['1x2motif + two-side\nbase-count beyond\n(12 params)']]

# create a plot for each productivity type
for (prod in c('productive', 'nonproductive')) {
    # subset data to productivity type
    subset = all_eval_results[productivity == prod]
    label_data_subset = label_data[productivity == prod]

    if (prod == 'productive'){
        bound = c(1.7, 2.71) 
    } else {
        bound = c(2, 2.73) 
    }

    # create plot
    require(ggrepel)
    plot = ggplot(subset) +
        geom_segment(x = 'full TCRB V-gene\ntraining dataset', y = motif_base_count_training_loss, xend = 'IGH testing\ndataset', yend = motif_base_count_training_loss, color = motif_base_count_color, linetype = 'dashed', size = 2) +
        geom_segment(x = 'TCRB training\ndataset J-genes\n(not included during\nmodel training)', y = motif_base_count_training_loss, xend = 'IGH testing\ndataset', yend = motif_base_count_training_loss, color = motif_base_count_color, linetype = 'dashed', size = 2) +
        geom_segment(x = 'TCRB training\ndataset V-genes\n(not included during\nmodel training)', y = motif_base_count_training_loss, xend = 'IGH testing\ndataset', yend = motif_base_count_training_loss, color = motif_base_count_color, linetype = 'dashed', size = 2) +
        geom_point(aes(y = loss, x = nice_loss_type, color = nice_model_type), size = 6)+
        geom_line(aes(y = loss, x = nice_loss_type, group = nice_model_type, color = nice_model_type), size = 4, alpha = 0.8)+
        geom_text_repel(data = label_data_subset, aes(y = loss, x = nice_loss_type, label = nice_model_type, color = nice_model_type), nudge_x = 0.2, fontface = "bold", size = 6, direction = 'y', hjust = 0, point.padding = 1, max.overlaps = Inf, lineheight = 0.8) +
        facet_wrap(~trim_type, nrow = 2, scales = 'free_x')+
        theme_cowplot(font_family = 'Arial') + 
        xlab(' ') +
        ylab('Expected per-sequence log loss\n')+
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.1) +
        theme(legend.position = 'none', text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), strip.text = element_text(size = 20), panel.spacing = unit(3, "lines"))  +
        scale_x_discrete(expand = expansion(add = c(0.2, 1.2)))+
        ylim(bound) +
        scale_color_manual(values = colors)

    # save plot
    ANNOTATION_TYPE <<- 'igor'
    path = get_manuscript_path()
    file_name = paste0(path, '/loss_compare_validation_', prod, '.pdf')
    ggsave(file_name, plot = plot, width = 16, height = 12, units = 'in', dpi = 750, device = cairo_pdf)
}

subset = subset[motif_type == MOTIF_TYPE]
cols = c('nice_model_type', 'trim_type', 'productivity', 'nice_loss_type', 'loss')
loss_data = subset[, ..cols]
colnames(loss_data) = c('model_type', 'trim_type', 'productivity', 'loss_type', 'log_loss')
fwrite(loss_data, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/validation_loss/loss.tsv'), sep = '\t')


