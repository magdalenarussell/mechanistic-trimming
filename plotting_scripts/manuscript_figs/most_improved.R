source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(plyr)
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

ANNOTATION_TYPE <<- 'igor' 

TRIM_TYPE <<- 'v_trim'

PRODUCTIVITY <<- 'nonproductive' 

MOTIF_TYPE <<- 'unbounded' 

NCPU <<- 2

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- 'all_subjects'

GENE_WEIGHT_TYPE <<- 'p_gene_marginal' 

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 0
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 0 

UPPER_TRIM_BOUND <<- as.numeric(14) 

MODEL_TYPE <<- 'two-side-base-count'

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(10)

RESIDUAL_COMPARE_FEATURE <<- NULL
TYPE <<- 'v_gene_family_loss'
source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/residual_comparison_functions.R')
source('analysis_scripts/pwm_profile_functions.R')
source('scripts/model_evaluation_functions.R')

genes = c('TRBV9', 'TRBV13')

# Read in dist data
predicted_trims1 = get_predicted_distribution_data() 
per_gene_resid1 = calculate_rmse_by_gene(predicted_trims1) 
per_gene_resid1$model = MODEL_TYPE
per_gene_trim_resid1 = calculate_rmse_by_gene_trim(predicted_trims1)
per_gene_trim_resid1$model = MODEL_TYPE
predicted_trims1$model = 'two-side-base-count model,\nall training data' 

MODEL_TYPE1 = MODEL_TYPE
LEFT_SIDE1 = LEFT_SIDE_TERMINAL_MELT_LENGTH
MODEL_TYPE <<- 'motif_two-side-base-count-beyond'
stopifnot(MODEL_TYPE %like% 'motif')
# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)


source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')
source('analysis_scripts/pwm_profile_functions.R')

# Read in dist data
predicted_trims2 = get_predicted_distribution_data() 
per_gene_resid2 = calculate_rmse_by_gene(predicted_trims2) 
per_gene_resid2$model = MODEL_TYPE
per_gene_trim_resid2 = calculate_rmse_by_gene_trim(predicted_trims2)
per_gene_trim_resid2$model = MODEL_TYPE
predicted_trims2$model = '+ 1x2 motif,\nall training data' 

#merge two predictions
together = rbind(per_gene_resid1, per_gene_resid2)
together_trim = rbind(per_gene_trim_resid1, per_gene_trim_resid2)
predictions = rbind(predicted_trims1, predicted_trims2, fill = TRUE)

together$model = factor(together$model, levels = c(MODEL_TYPE1, MODEL_TYPE))

wide = together[, -c('p_gene')] %>% pivot_wider(names_from = 'model', values_from = 'rmse') %>% as.data.table()
slopes = wide[,(get(MODEL_TYPE) - get(MODEL_TYPE1)), by = gene]
setnames(slopes, 'V1', 'slope')
together = merge(together, slopes)

together = together[order(abs(slope))]

outlier_count = ceiling(nrow(together)*0.1)
if ((outlier_count %% 2) != 0){
    outlier_count = outlier_count + 1
}
outliers = together[(nrow(together)-outlier_count + 1):nrow(together)]
no_outliers = together[1:(nrow(together)-outlier_count)]

cutoff = mean(c(min(no_outliers$slope), max(outliers$slope)))

wide$diff = wide[['motif_two-side-base-count-beyond']] - wide[['two-side-base-count']]
plot_data = wide
plot = ggplot(plot_data) +
    geom_histogram(aes(x = diff), binwidth= 0.01) +
    geom_vline(xintercept = cutoff, size = 5, color = 'gray') +
    geom_text(x = -0.15, y = 6, label = 'improved\ngenes', color = 'gray', size = 12, lineheight = 0.8, check_overlap = TRUE)+
    ylab('Gene count\n') +
    xlab('\nPer-gene RMSE difference') +
    theme_cowplot(font_family = 'Arial') +
    theme(text = element_text(size = 30), axis.text.x=element_text(size = 30), axis.text.y = element_text(size = 30), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

path = get_manuscript_path()
file_name = paste0(path, '/most_improved.pdf')
ggsave(file_name, plot = plot, width = 15, height = 9, units = 'in', dpi = 750, device = cairo_pdf)

# re-fit_model without outliers
motif_data = aggregate_all_subject_data()
new_motif_data = motif_data[gene %in% no_outliers$gene]

# plot predicted distributions on held out (outlier genes)
model = fit_model(new_motif_data)
outlier_motif_data = motif_data[gene %in% outliers$gene]
outlier_motif_data = process_data_for_model_fit(outlier_motif_data)
outlier_motif_data$predicted_prob = temp_predict(model, newdata = outlier_motif_data)
outlier_motif_data[, empirical_prob := count/sum(count), by = .(subject, gene)]
outlier_motif_data$model = paste0('+ 1x2motif,\nremoving improved genes')

predictions = rbind(predictions, outlier_motif_data, fill = TRUE)

# re-fit_model without outliers
gene_families = get_gene_families(cluster_count = 5, combine_by_terminal = FALSE, full_sequence = TRUE, align = TRUE)$cluster_data
outlier_clusters = unique(gene_families[gene %in% unique(outliers$gene)]$clusters_grouped)
outlier_cluster_genes = unique(gene_families[clusters_grouped %in% outlier_clusters]$gene)
no_outliers = together[!(gene %in% outliers) & !(gene %in% outlier_cluster_genes)]
new_motif_data = motif_data[gene %in% no_outliers$gene]

# plot predicted distributions on held out (outlier genes)
model = fit_model(new_motif_data)
outlier_motif_data = motif_data[gene %in% outliers$gene]
outlier_motif_data = process_data_for_model_fit(outlier_motif_data)
outlier_motif_data$predicted_prob = temp_predict(model, newdata = outlier_motif_data)
outlier_motif_data[, empirical_prob := count/sum(count), by = .(subject, gene)]

outlier_motif_data$model = paste0('+ 1x2motif,\nremoving improved genes + similar')

predictions = rbind(predictions, outlier_motif_data, fill = TRUE)

model_types = c('two-side-base-count model,\nall training data', '+ 1x2 motif,\nall training data', '+ 1x2motif,\nremoving improved genes', '+ 1x2motif,\nremoving improved genes + similar')

predictions$model = factor(predictions$model, levels = model_types) 

for (gene_name in genes){
    important_cols = c('trim_length', 'predicted_prob', 'gene', 'model')
    predicted_data = unique(predictions[gene == gene_name, ..important_cols])
    predicted_data$model = factor(predicted_data$model, levels = model_types)
    empirical_data = predictions[gene == gene_name][order(subject, trim_length)]
    empirical_data$model = factor(empirical_data$model, levels = model_types)
    title = paste0(gene_name)
    # get gene sequence
    gene_seq = get_gene_sequence(gene_name, max(predictions$trim_length))
    gene_seq_with_positions = get_plot_positions_for_gene_sequence(gene_seq)
    
    max_prob = max(max(empirical_data$empirical_prob), max(predicted_data$predicted_prob))

    labels = data.table(model = model_types, yvar = max_prob, leftx = -2.1, rightx = UPPER_TRIM_BOUND + 0.1) 
    temp_plot = ggplot() +
        geom_line(data = empirical_data, aes(x = trim_length, y = empirical_prob, group = subject), size = 1, alpha = 0.5, color = 'grey') +
        geom_line(data = predicted_data, aes(x = trim_length, y = predicted_prob), size = 2, alpha = 0.7, color = 'blue') +
        facet_grid(cols = vars(factor(model, levels = model_types))) +
        geom_vline(xintercept = 0, color = 'black', size = 3) +
        geom_text(data = gene_seq_with_positions, y = max_prob, aes(x = position, label = base), size = 7) +
        geom_text(data = labels, aes(y = yvar, x = leftx), label = '3\'- ', size = 6) +
        geom_text(data = labels, aes(y = yvar, x = rightx), label = ' -5\'', size = 6) +
        ggtitle(title) +
        xlab('Number of trimmed nucleotides') +
        ylab('Probability\n') +
        theme_cowplot(font_family = 'Arial') + 
        theme(legend.position = "none", text = element_text(size = 30), axis.text.x=element_text(size = 25), axis.text.y = element_text(size = 25), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        ylim(c(0, max_prob + 0.02))
    assign(gene_name, temp_plot)
} 

all = align_plots(TRBV9, TRBV13, plot, align = 'v', axis = 'l')

first_row = plot_grid(all[[3]], NULL, nrow = 1, rel_widths = c(1, 0))

all_tog = plot_grid(first_row, NULL, all[[1]], NULL, all[[2]], NULL, nrow = 6, rel_heights = c(0.8, 0.05, 0.5, 0.05, 0.5, 0.05), labels = c('A', '', 'B', '', 'C', ''), label_size = 35)

path = get_manuscript_path()
file_name = paste0(path, '/motif_exploration.pdf')
ggsave(file_name, plot = all_tog, width = 27, height = 23, units = 'in', dpi = 750, device = cairo_pdf)


