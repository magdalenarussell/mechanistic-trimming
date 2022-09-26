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
LOSS_GENE_WEIGHT <<- 'uniform' 
stopifnot(LOSS_GENE_WEIGHT %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform', 'p_gene_marginal_all_seqs', 'p_gene_given_subject_all_seqs'))

source('scripts/model_evaluation_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/model_evaluation_functions.R')

all_eval_results = data.table()
# annotation_types = c('validation_data_alpha', 'validation_data_beta', 'validation_data_gamma')
annotation_types = c('validation_data_alpha', 'validation_data_beta', 'igor')

trim_types = c('j_trim', 'v_trim')
prods = c('productive', 'nonproductive')

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
            setnames(temp_eval_results, type, 'loss')
            temp_eval_results$loss_type = type
            temp_eval_results$trim_type = trim_type
            temp_eval_results$productivity = prod
            all_eval_results = rbind(all_eval_results, temp_eval_results, fill = TRUE)
        }
    }
}

# get model types
orig_model_types = c("motif_two-side-base-count-beyond", "motif", 'null')
new_model_types = c('1x2motif + two-side\nbase-count beyond\n(12 params)', '2x4motif (18 params)', 'null (0 params)')

all_eval_results$model_type = mapvalues(all_eval_results$model_type, from = orig_model_types, to = new_model_types) 
all_eval_results[loss_type == 'igor', loss_type := paste0(loss_type, '_', trim_type, '_', productivity)]
neat_names = make_model_names_neat(new_model_types)
colors = set_color_palette(c(neat_names, '2x4motif (18 params)'), with_params = TRUE)

loss_types = unique(all_eval_results$loss_type)

nice_loss_types = c('TCRA testing\ndataset', 'TCRB testing\ndataset', 'TCRB training\ndataset J-genes\n(not included during\nmodel training)', 'TCRB training\ndataset J-genes\n(not included during\nmodel training)', 'TCRB training\ndataset V-genes\n(not included during\nmodel training)', 'full TCRB V-gene\ntraining dataset')
loss_order = c('TCRB training\ndataset J-genes\n(not included during\nmodel training)', 'full TCRB V-gene\ntraining dataset', 'TCRB training\ndataset V-genes\n(not included during\nmodel training)', 'TCRB testing\ndataset', 'TCRA testing\ndataset')
all_eval_results$nice_loss_type = mapvalues(all_eval_results$loss_type, from = loss_types, to=nice_loss_types)
all_eval_results[nice_loss_type == 'full TCRB V-gene\ntraining dataset', trim_type := 'v_trim']
all_eval_results$trim_type = mapvalues(all_eval_results$trim_type, from = unique(all_eval_results$trim_type), to=c('J-gene trimming', 'V-gene trimming'))
all_eval_results$trim_type = factor(all_eval_results$trim_type, levels = c('V-gene trimming', 'J-gene trimming'))

# bound = c(2.04, 3.01)
bound = c(1.9, 2.7)

# reformat model names to be nice
all_eval_results = all_eval_results[motif_type == MOTIF_TYPE]
nice_names = make_model_names_neat(unique(all_eval_results$model_type)) 
all_eval_results$nice_model_type = mapvalues(all_eval_results$model_type, from = unique(all_eval_results$model_type), to = nice_names)

ordered_losses = loss_order
all_eval_results$nice_loss_type = factor(all_eval_results$nice_loss_type, levels = ordered_losses)
last_loss = ordered_losses[length(ordered_losses)]
label_data = all_eval_results[nice_loss_type == last_loss] 

motif_base_count_training_loss = all_eval_results[nice_loss_type == 'full TCRB V-gene\ntraining dataset' & nice_model_type == '1x2motif + two-side\nbase-count beyond\n(12 params)']$loss
motif_base_count_color = colors[['1x2motif + two-side\nbase-count beyond\n(12 params)']]

for (prod in c('productive', 'nonproductive')) {
    subset = all_eval_results[productivity == prod]
    label_data_subset = label_data[productivity == prod]
    # create plot
    require(ggrepel)
    plot = ggplot(subset) +
        geom_segment(x = 'full TCRB V-gene\ntraining dataset', y = motif_base_count_training_loss, xend = 'TCRA testing\ndataset', yend = motif_base_count_training_loss, color = motif_base_count_color, linetype = 'dashed', size = 2) +
        geom_segment(x = 'TCRB training\ndataset J-genes\n(not included during\nmodel training)', y = motif_base_count_training_loss, xend = 'TCRA testing\ndataset', yend = motif_base_count_training_loss, color = motif_base_count_color, linetype = 'dashed', size = 2) +
        geom_segment(x = 'TCRB training\ndataset V-genes\n(not included during\nmodel training)', y = motif_base_count_training_loss, xend = 'TCRA testing\ndataset', yend = motif_base_count_training_loss, color = motif_base_count_color, linetype = 'dashed', size = 2) +
        geom_point(aes(y = loss, x = nice_loss_type, color = nice_model_type), size = 12)+
        geom_line(aes(y = loss, x = nice_loss_type, group = nice_model_type, color = nice_model_type), size = 7)+
        geom_text_repel(data = label_data_subset, aes(y = loss, x = nice_loss_type, label = nice_model_type, color = nice_model_type), nudge_x = 0.2, fontface = "bold", size = 12, direction = 'y', hjust = 0, point.padding = 1, max.overlaps = Inf, lineheight = 0.8) +
        facet_wrap(~trim_type, nrow = 2, scales = 'free_x')+
        theme_cowplot(font_family = 'Arial') + 
        xlab(' ') +
        ylab('Expected per-sequence log loss (uniform gene usage)\n')+ 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.1) +
        theme(legend.position = 'none', text = element_text(size = 38), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 38), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), strip.text = element_text(size = 38), panel.spacing = unit(3, "lines"))  +
        scale_x_discrete(expand = expansion(add = c(0.2, 1.8)))+
        ylim(bound) +
        scale_color_manual(values = colors)


    ANNOTATION_TYPE <<- 'igor'
    path = get_manuscript_path()
    file_name = paste0(path, '/loss_compare_validation_', prod, '_uniform_gene_loss.pdf')
    if (prod == 'nonproductive'){
        w =20 
    } else {
        w = 20 
    }
    ggsave(file_name, plot = plot, width = w, height = 25, units = 'in', dpi = 750, device = cairo_pdf)
}


