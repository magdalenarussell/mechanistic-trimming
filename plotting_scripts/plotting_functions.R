source(paste0('plotting_scripts/model_group_functions/', MODEL_GROUP, '.R'))
source(paste0(PROJECT_PATH, '/plotting_scripts/plot_paths.R'))

get_gene_sequence <- function(gene_name, gene_seq_length, pnuc_count = 2){
    whole_nucseq = get_oriented_whole_nucseqs()
    temp_data = whole_nucseq[toupper(substring(gene, 4, 4)) == toupper(substring(GENE_NAME, 1,1))]
    setnames(whole_nucseq, 'gene', GENE_NAME)
    colnames(temp_data) = c(GENE_NAME, 'sequence')
    gene_groups = get_common_genes_from_seqs(temp_data)
    together = unique(merge(temp_data, gene_groups)[, c('gene', 'sequence')]) 
    gene = together[gene == gene_name][1]

    # get sequence
    whole_gene_seq = DNAString(gene$sequence)
    possible_pnucs_5_to_3 = substring(reverseComplement(whole_gene_seq),1, pnuc_count)
    whole_gene_with_pnucs = c(unlist(whole_gene_seq), unlist(possible_pnucs_5_to_3))
    subset = substring(whole_gene_with_pnucs, nchar(whole_gene_with_pnucs) - (gene_seq_length + pnuc_count-1), nchar(whole_gene_with_pnucs))
    return(subset)
}

get_plot_positions_for_gene_sequence <- function(gene_sequence, pnuc_count = 2){
    end_position = -1 * pnuc_count + 0.5 
    start_position = nchar(gene_sequence) - pnuc_count - 0.5
    positions = seq(start_position, end_position, by = -1)
    together = data.table(base = str_split(unlist(gene_sequence), '')[[1]], position = positions)
    return(together)
}

plot_predicted_trimming_dists_single_group <- function(data, gene_name, file_name){
    important_cols = c('trim_length', 'predicted_prob', 'gene')
    predicted_data = data[gene == gene_name, ..important_cols]
    empirical_data = data[gene == gene_name][order(subject, trim_length)]
    total_tcrs = paste(data[gene == gene_name, sum(count)])
    title = paste0(gene_name, ', Total TCR count = ', total_tcrs)
    # get gene sequence
    gene_seq = get_gene_sequence(gene_name, max(data$trim_length))
    gene_seq_with_positions = get_plot_positions_for_gene_sequence(gene_seq)
    
    max_prob = max(max(empirical_data$empirical_prob), max(predicted_data$predicted_prob))

    plot = ggplot() +
        geom_line(data = empirical_data, aes(x = trim_length, y = empirical_prob, group = subject), size = 2, alpha = 0.7, color = 'grey') +
        geom_line(data = predicted_data, aes(x = trim_length, y = predicted_prob), size = 3, alpha = 0.7, color = 'blue') +
        geom_vline(xintercept = 0, color = 'black', size = 4) +
        geom_text(data = gene_seq_with_positions, y = max_prob, aes(x = position, label = base), size = 8) +
        geom_text(y = max_prob, aes(x = -2.1), label = '3\' -', size = 10) +
        geom_text(y = max_prob, aes(x = UPPER_TRIM_BOUND + 0.1), label = '- 5\'', size = 10) +
        ggtitle(title) +
        xlab('Number of trimmed nucleotides') +
        ylab('Probability') +
        theme_cowplot(font_family = 'Arial') + 
        theme(legend.position = "none", text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 

    ggsave(file_name, plot = plot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
    print(paste0('finished plot for ', gene_name))
    return(plot)
}

map_positions_to_values <- function(positions){
    values = c()
    for (position in positions){
        if (substring(position, 1, 7) == "motif_5"){
            val = -1 * as.numeric(substring(position, nchar(position), nchar(position)))
        } else if (substring(position, 1, 7) == "motif_3"){
            val = as.numeric(substring(position, nchar(position), nchar(position)))
        }
        values = c(values, val)
    }
    together = data.table(positions = positions, values = values)
    return(together)
}

plot_melting_coefficient_heatmap_single_group <- function(model_coef_matrix, file_name, with_values = FALSE, limits = NULL, write_plot = TRUE){
    model_coef_matrix = model_coef_matrix[parameter %like% 'terminal_melting']
    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$coefficient/log(10)
    # order variables
    model_coef_matrix$parameter = factor(model_coef_matrix$parameter, levels = c('left_terminal_melting', 'right_terminal_melting'))

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    plot = ggplot(model_coef_matrix, aes(x=parameter, y=base, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Melting temperature') +
        ylab('') +
        geom_vline(xintercept = 1 + 0.5, size = 3, color = 'black') +
        theme(text = element_text(size = 15), axis.text.y = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
        guides(fill = guide_colourbar(barheight = 8)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits) +
        annotate("text", x = 0.55, y = 0.45, label = "5\'", size = 6) +  
        annotate("text", x = 2.55 , y = 0.45, label = "3\'", size = 6) 
    
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = model_coef_matrix, aes(x = parameter, y = base, label = round(log_10_pdel, 2)))
    }

    if (isTRUE(write_plot)){
        ggsave(file_name, plot = plot, width = 7.5, height = 3, units = 'in', dpi = 750, device = cairo_pdf)
    } else {
        return(plot)
    }
}


plot_distance_coefficient_heatmap_single_group <- function(model_coef_matrix, file_name, with_values = FALSE, limits = NULL, write_plot = TRUE){
    model_coef_matrix = model_coef_matrix[parameter %like% 'trim_length']
    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$coefficient/log(10)
    # order variables
    model_coef_matrix$parameter = factor(model_coef_matrix$parameter, levels = paste0('trim_length_', seq(UPPER_TRIM_BOUND, LOWER_TRIM_BOUND)))
    total_dist = UPPER_TRIM_BOUND - LOWER_TRIM_BOUND

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    plot = ggplot(model_coef_matrix, aes(x=parameter, y=base, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Distance') +
        ylab('') +
        theme(text = element_text(size = 35), axis.text.x = element_text(size = 22), axis.text.y = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
        guides(fill = guide_colourbar(barheight =12)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits) 
    
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = model_coef_matrix, aes(x = parameter, y = base, label = round(log_10_pdel, 2)), size = 10)
    }

    if (isTRUE(write_plot)){
        ggsave(file_name, plot = plot, width = 3*total_dist, height = 4, units = 'in', dpi = 750, device = cairo_pdf)
    } else {
        return(plot)
    }
}

plot_model_coefficient_heatmap_single_group <- function(model_coef_matrix, file_name, with_values = FALSE, limits = NULL, write_plot = TRUE){
    model_coef_matrix = model_coef_matrix[parameter %like% 'motif']

    position_values = map_positions_to_values(unique(model_coef_matrix$parameter))
    together = merge(model_coef_matrix, position_values, by.x = 'parameter', by.y = 'positions')

    # convert to log_10
    together$log_10_pdel = together$coefficient/log(10)
    # order variables
    together$base = factor(together$base, levels = c('T', 'G', 'C', 'A'))
    together$values = factor(together$values)

    if (is.null(limits)){
        max_val = max(abs(together$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    motif_length = RIGHT_NUC_MOTIF_COUNT + LEFT_NUC_MOTIF_COUNT

    plot = ggplot(together, aes(x=values, y=base, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Position') +
        ylab ('Base') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        geom_vline(xintercept = LEFT_NUC_MOTIF_COUNT + 0.5, size = 3, color = 'black') +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits) +
        annotate("text", x = 0.35, y = 0.25, label = "5\'", size = 8) +  
        annotate("text", x = motif_length + 0.65, y = 0.25, label = "3\'", size = 8) +  
        coord_cartesian(ylim = c(1, 4), clip = "off")
    
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = together, aes(x = values, y = base, label = round(log_10_pdel, 2)))
    }

    if (isTRUE(write_plot)){
        ggsave(file_name, plot = plot, width = motif_length + 4, height = 4, units = 'in', dpi = 750, device = cairo_pdf)
    } else {
        return(plot)
    }
}

plot_background_base_composition_heatmap_single_group <- function(background_matrix, file_name, with_values = FALSE, write_plot = TRUE){
    data = background_matrix %>%
        pivot_longer(!c(base), names_to = 'position', values_to = 'coef')

    position_values = map_positions_to_values(unique(data$position))
    together = merge(data, position_values, by.x = 'position', by.y = 'positions')

    # order variables
    together$base = factor(together$base, levels = c('T', 'G', 'C', 'A'))
    together$values = factor(together$values)

    motif_length = RIGHT_NUC_MOTIF_COUNT + LEFT_NUC_MOTIF_COUNT
    
    together$log_fold_freq = log2(together$coef) - log2(0.25)

    max = round(max(abs(together$log_fold_freq)), 1) + 0.1 
    plot = ggplot(together, aes(x=values, y=base, fill=log_fold_freq)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Position') +
        ylab ('Base') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        geom_vline(xintercept = LEFT_NUC_MOTIF_COUNT + 0.5, size = 3, color = 'black') +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log2 fold change in\n background base frequency', limits = c(-1*max, max)) +
        annotate("text", x = 0.35, y = 0.25, label = "5\'", size = 8) +  
        annotate("text", x = motif_length + 0.65, y = 0.25, label = "3\'", size = 8) +  
        coord_cartesian(ylim = c(1, 4), clip = "off")
    
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = together, aes(x = values, y = base, label = round(log_fold_freq, 2)))
    }

    if (isTRUE(write_plot)){
        ggsave(file_name, plot = plot, width = 9, height = 4.1, units = 'in', dpi = 750, device = cairo_pdf)
    } else {
        return(plot)
    }
}

plot_model_residual_boxplot_single_group <- function(data, gene_name, file_name){
    data$residual = data$empirical_prob - data$predicted_prob
    data_subset = data[gene == gene_name]
    total_tcrs = paste(data[gene == gene_name, sum(count)])
    title = paste0(gene_name, ', Total TCR count = ', total_tcrs)
 
    plot = ggplot(data = data_subset, aes(x = as.factor(trim_length), y = residual)) +
        geom_boxplot(width=0.8, color="black", size = 1.2) +
        geom_hline(yintercept = 0, color = 'gray', size = 3) +
        xlab('Number of trimmed nucleotides') +
        ylab('Observed prob - Predicted prob') +
        ggtitle(title) +
        theme_cowplot(font_family = 'Arial') + 
        theme(legend.position = "none", text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 
    
    ggsave(file_name, plot = plot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
    print(paste0('finished plot for ', gene_name))
    return(plot) 
}

# let model_coef_matrix contain coefficient data for individual subjects
plot_model_coefficient_variations <- function(model_coef_matrix){
    all_subject_data = get_all_subject_model_coefficient_data()
    combined = rbind(model_coef_matrix, all_subject_data)
    data = combined %>%
        pivot_longer(!c(base, model_group), names_to = 'position', values_to = 'coef')

    position_values = map_positions_to_values(unique(data$position))
    together = merge(data, position_values, by.x = 'position', by.y = 'positions')

    # convert to log_10
    together$log_10_pdel = together$coef/log(10)
    # order variables
    together$base = factor(together$base, levels = c('A', 'C', 'G', 'T'))
    together$values = factor(together$values)
    together = data.table(together)

    plot = ggplot(together[model_group !='all_subjects'], aes(x=log_10_pdel)) +
        facet_grid(rows = vars(base), cols = vars(values)) +
        geom_histogram(aes(y = stat(width*density), group = values), position = 'identity')+
        geom_vline(data = together[model_group =='all_subjects'], aes(xintercept=log_10_pdel), color = '#1b9e77', linetype = 'F1', size=1) +
        theme_cowplot(font_family = 'Arial') + 
        xlab('log10(probability of deletion)') +
        ylab ('Proportion of subjects') +
        theme(text = element_text(size = 20)) +
        # geom_vline(xintercept = RIGHT_NUC_MOTIF_COUNT + 0.5, size = 3, color = 'black') +
        background_grid(major = 'xy') +
        panel_border(color = 'gray60', size = 1.5)

    path = get_coef_variations_file_path()
    file = 'coefficient_variations.pdf'

    file_name = file.path(path, file)
    ggsave(file_name, plot = plot, width = 15, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_all_model_residuals_hist <- function(data, file_name){
    data$residual = data$empirical_prob - data$predicted_prob
    mean_data = data[, mean(abs(residual)), by = .(gene)]
    setnames(mean_data, 'V1', 'mean_residual')
    overall_mean_data = data[, mean(abs(residual))]

    plot = ggplot() +
        geom_histogram(data = mean_data, mapping = aes(x = mean_residual), alpha = 0.75) +
        geom_vline(xintercept = 0, color = 'gray', size = 3) +
        geom_vline(xintercept = overall_mean_data, color = 'blue', size = 2, alpha = 0.8) +
        xlab('abs(Observed prob - Predicted prob)') +
        theme_cowplot(font_family = 'Arial') + 
        theme(legend.position = "none", text = element_text(size = 30), axis.text = element_text(size = 20), axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 
    ggsave(file_name, plot = plot, width = 14, height = 8, units = 'in', dpi = 750, device = cairo_pdf)
    return(plot) 
}

plot_all_model_residuals_plot <- function(data, file_name, color_gene_list = NULL){
    data$residual = data$empirical_prob - data$predicted_prob
    mean_data = data[, mean(residual), by = .(gene, trim_length)]
    setnames(mean_data, 'V1', 'mean_residual')
    overall_mean_data = data[, mean(residual), by = .(trim_length)]
    setnames(overall_mean_data, 'V1', 'overall_mean_residual')

    # var_data = data[, var(residual), by = .(trim_length)]
    # setnames(var_data, 'V1', 'var')

    if (!is.null(color_gene_list)){
        mean_data[, outlier := FALSE]
        mean_data[gene %in% color_gene_list, outlier := TRUE]
        plot = ggplot() +
            geom_line(data = mean_data, aes(x = trim_length, y = mean_residual, group = gene, color = outlier), size = 2, alpha = 0.5) 
    } else {
        plot = ggplot() +
            geom_line(data = mean_data, aes(x = trim_length, y = mean_residual, group = gene), size = 2, alpha = 0.5) 
    }
    
    plot = plot +
        geom_hline(yintercept = 0, color = 'gray', size = 3) +
        geom_line(data = overall_mean_data, aes(x = trim_length, y = overall_mean_residual), color = 'blue', size = 2, alpha = 0.8) +
        ylab('Observed prob - Predicted prob\n') +
        theme_cowplot(font_family = 'Arial') + 
        theme(legend.position = "none", text = element_text(size = 30), axis.text = element_text(size = 20), axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        ylim(-0.5, 1)
    ggsave(file_name, plot = plot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
    return(plot) 
}

plot_model_evaluation_heatmap <- function(eval_data, type, with_values = FALSE, model_type, terminal_melting_5_end_length_filter, limits = NULL){
    stopifnot(length(model_type) == 1)
    path = get_model_eval_file_path(type)
    if (model_type %like% 'two_side_terminal_melting') {
        model_type_name = paste0(model_type, '_', terminal_melting_5_end_length_filter)
    } else {
        model_type_name = model_type
    }
    file_name = paste0(path, '/model_evaluation_heatmap_model_type_filter_', model_type_name, '.pdf')
    eval_data = process_model_evaluation_file(eval_data, model_types_neat = model_type, terminal_melting_5_end_length_filter = terminal_melting_5_end_length_filter)
    eval_data = eval_data[motif_type == MOTIF_TYPE]

    if (type == 'expected_log_loss'){
        fill_lab = 'Expected conditional log loss'
    } else if (type == 'log_loss'){
        fill_lab = 'Conditional log loss\n(training dataset)'
    } else if (type == 'aic'){
        fill_lab = 'AIC' 
    } else if (type == 'raw_loss'){
        fill_lab = 'Raw count log loss\n(training dataset)' 
    } else if (type == 'old_loss_cv'){
        fill_lab = 'Raw count conditional log loss\n(held-out dataset)' 
    }

    setnames(eval_data, type, 'loss')

    plot = ggplot(eval_data) +
        geom_tile(aes(x = motif_length_5_end, y = motif_length_3_end, fill = loss)) +
        theme_cowplot(font_family = 'Arial') + 
        xlab('5\' motif length') +
        ylab('3\' motif length') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        guides(fill = guide_colourbar(barheight = 14, limits = limits)) +
        scale_fill_viridis_c(name = fill_lab) +
        scale_x_continuous(breaks = seq(0, 6)) +
        scale_y_continuous(breaks = seq(0, 6))
    
    if (!is.null(limits)){
        plot = plot + scale_fill_viridis_c(name = fill_lab, limits = limits) 
    } else {
        plot = plot + scale_fill_viridis_c(name = fill_lab) 
    }

    if (with_values == TRUE){
        plot = plot +
            geom_text(data = eval, aes(x = motif_length_5_end, y = motif_length_3_end, label = round(loss, 0)))
    }

    ggsave(file_name, plot = plot, width = 11, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_model_evaluation_scatter_coef_count <- function(eval_data, type, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, label = FALSE) {
    # process evaluation file
    eval_data = process_model_evaluation_file(eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter)
    eval_data = eval_data[motif_type == MOTIF_TYPE]

    if (type == 'expected_log_loss'){
        ylab = 'Expected conditional log loss'
    } else if (type == 'log_loss'){
        ylab = 'Conditional log loss (training dataset)'
    } else if (type == 'aic'){
        ylab = 'AIC' 
    } else if (type == 'raw_loss'){
        ylab = 'Raw count log loss (training dataset)' 
    } else if (type == 'old_loss_cv'){
        ylab = 'Raw count conditional log loss\n(held-out dataset)' 
    }

    plot = ggplot(eval_data) +
        geom_point(aes(y = get(type), x = model_parameter_count, color = model_type), size = 5)+
        theme_cowplot(font_family = 'Arial') + 
        xlab('Total number of terms') +
        ylab(ylab) +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 

    path = get_model_eval_file_path(type)
    file_name = paste0(path, '/neat_', type, '_term_count_scatter_', left_motif_size_filter, '_', right_motif_size_filter, '_motif')

    if (isTRUE(label)){
        plot = plot +
            geom_text_repel(aes(y = get(type), x = model_parameter_count, color = model_type, label = model_type), size = 4) + 
            scale_x_continuous(breaks = seq(0, 60, 2), limits = c(0, 60)) +
            theme(legend.position = 'none', text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) 
        file_name = paste0(file_name, '.pdf')
    } else {
        plot = plot +
            theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) 
        file_name = paste0(file_name, '_no_label.pdf')
    }
    ggsave(file_name, plot = plot, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_model_evaluation_compare <- function(eval_data1, eval_data2, type1, type2, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, label = FALSE) {
    # merge two types
    cols = colnames(eval_data1)
    cols = cols[!(cols %in% c('held_out_gene_fraction', 'sample_repetitions'))]
    eval_data = merge(eval_data1, eval_data2, by = cols[!(cols == type1)])

    # process evaluation file
    eval_data = process_model_evaluation_file(eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter)
    eval_data = eval_data[motif_type == MOTIF_TYPE]

    if (type1 == 'expected_log_loss'){
        xlab = 'Expected conditional log loss'
    } else if (type1 == 'log_loss'){
        xlab = 'Conditional log loss (training dataset)'
    } else if (type1 == 'aic'){
        xlab = 'AIC' 
    } else if (type1 == 'raw_loss'){
        xlab = 'Raw count log loss (training dataset)' 
    } else if (type1 == 'old_loss_cv'){
        xlab = 'Raw count conditional log loss\n(held-out dataset)' 
    }

    if (type2 == 'expected_log_loss'){
        ylab = 'Expected conditional log loss'
    } else if (type2 == 'log_loss'){
        ylab = 'Conditional log loss (training dataset)'
    } else if (type2 == 'aic'){
        ylab = 'AIC' 
    } else if (type2 == 'raw_loss'){
        ylab = 'Raw count log loss (training dataset)' 
    } else if (type2 == 'old_loss_cv'){
        ylab = 'Raw count conditional log loss\n(held-out dataset)' 
    }
    
    plot = ggplot(eval_data) +
        geom_abline(intercept = 0, slope = 1, size = 3, color = 'gray60', linetype = 'dashed') +
        geom_point(aes(y = get(type2), x = get(type1), color = model_type), size = 5)+
        theme_cowplot(font_family = 'Arial') + 
        xlab(xlab) +
        ylab(ylab) +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 

    path = get_model_eval_file_path('compare')
    file_name = paste0(path, '/neat_', type1, '_', type2, '_eval_compare_', left_motif_size_filter, '_', right_motif_size_filter, '_motif')

    if (isTRUE(label)){
        plot = plot +
            geom_text_repel(aes(y = get(type2), x = get(type1), color = model_type, label = model_type), size = 4) + 
            scale_x_continuous(breaks = seq(0, 60, 2), limits = c(0, 60)) +
            theme(legend.position = 'none', text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) 
        file_name = paste0(file_name, '.pdf')
    } else {
        plot = plot +
            theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) 
        file_name = paste0(file_name, '_no_label.pdf')
    }
    ggsave(file_name, plot = plot, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_model_evaluation_compare_bounds <- function(eval_data, upper_bound1, upper_bound2, type, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, label = FALSE, limits = NULL) {
    # process evaluation file
    eval_data1 = process_model_evaluation_file(eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, upper_trim_bound = upper_bound1)
    eval_data2 = process_model_evaluation_file(eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, upper_trim_bound = upper_bound2)

    cols = colnames(eval_data1)
    together = merge(eval_data1, eval_data2, by = cols[!(cols %in% c('upper_bound', 'model_parameter_count', type))])
    together = together[motif_type == MOTIF_TYPE]

    if (type == 'expected_log_loss'){
        lab = 'Expected conditional log loss'
    } else if (type == 'log_loss'){
        lab = 'Conditional log loss (training dataset)'
    } else if (type == 'aic'){
        lab = 'AIC' 
    } else if (type == 'raw_loss'){
        lab = 'Raw count log loss (training dataset)' 
    } else if (type == 'old_loss_cv'){
        lab = 'Raw count conditional log loss\n(held-out dataset)' 
    }

    plot = ggplot(together) +
        geom_abline(intercept = 0, slope = 1, size = 3, color = 'gray60', linetype = 'dashed') +
        geom_point(aes(y = get(paste0(type, '.y')), x = get(paste0(type, '.x')), color = model_type), size = 5)+
        theme_cowplot(font_family = 'Arial') + 
        xlab(paste0(lab, '\n trims bounded at ', upper_bound1)) +
        ylab(paste0(lab, '\n trims bounded at ', upper_bound2)) +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 
    
    if (!is.null(limits)){
        plot = plot +
            scale_x_continuous(limits = limits) +
            scale_y_continuous(limits = limits)
    }

    path = get_model_eval_file_path('compare_bounds')
    file_name = paste0(path, '/neat_', type, '_eval_compare_bounds_', upper_bound1, '_and_', upper_bound2, '_motifs_', left_motif_size_filter, '_', right_motif_size_filter)

    if (isTRUE(label)){
        plot = plot +
            geom_text_repel(aes(y = get(paste0(type, '.y')), x = get(paste0(type, '.x')), color = model_type, label = model_type), size = 4) + 
            scale_x_continuous(breaks = seq(0, 60, 2), limits = c(0, 60)) +
            theme(legend.position = 'none', text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) 
        file_name = paste0(file_name, '.pdf')
    } else {
        plot = plot +
            theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) 
        file_name = paste0(file_name, '_no_label.pdf')
    }
    ggsave(file_name, plot = plot, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
}

get_base_composition_counts <- function(motif_data){
    positions = get_positions()
    processed = data.table(base = c('A', 'T', 'C', 'G'))
    for (pos in positions){
        temp = motif_data[, .N, by = get(pos)]
        temp[, prop := N/sum(N)]
        colnames(temp) = c('base', 'count', pos)
        processed = merge(processed, temp[, -c('count')], by = 'base')
    }
    return(processed)
}

plot_coefficient_variation_across_individuals <- function(all_individual_coefficients){
    if (!('base' %in% colnames(all_individual_coefficients))){
            all_individual_coefficients$base = '' 
    }
    all_individual_coefficients[!(base ==  ''), parameter := paste0(parameter, ', base ', base)]

    var_data = all_individual_coefficients[, var(coefficient), by = parameter]
    setnames(var_data, 'V1', 'coef_var')
    var_data$parameter = factor(var_data$parameter, levels = var_data[order(coef_var)]$parameter)

    plot = ggplot(var_data) +
        geom_point(aes(x = coef_var, y = parameter), size = 5) +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Coefficient variance') +
        ylab('Parameter') +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        theme(legend.position = 'none', text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) 
    
    file_path = get_individual_comparison_file_path()
    file_name = paste0(file_path, '/coefficient_variance.pdf')
    ggsave(file_name, plot = plot, width = 14, height = 16, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_coefficient_by_snp <- function(coef_snp_data, snpID, parameter_group = NULL){
    require(ggpubr)
    file_path = get_individual_comparison_file_path()
    if (!is.null(parameter_group)){
        coef_snp_data = coef_snp_data[parameter %like% parameter_group]
        file_name = paste0(file_path, '/', parameter_group, '_parameters_by_snp', snpID, '.pdf')
    } else {
        file_name = paste0(file_path, '/parameters_by_snp', snpID, '.pdf')
    }
    if (parameter_group == 'trim_length'){
        coef_snp_data$parameter = factor(coef_snp_data$parameter, levels = paste0('trim_length_', seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)))
    }
    subset = coef_snp_data[snp == snpID & !is.na(genotype)]
    subset$genotype = as.character(subset$genotype)
    subset$log_10_coef = subset$coefficient/log(10)
    comparisons = list(c("0", "1"), c("1", "2"), c("0", "2")) 
    row_count = ceiling(length(unique(subset$parameter))/6)
    col_count = min(6, length(unique(subset$parameter)))
    plot = ggplot(subset, aes(x = genotype, y = log_10_coef)) +
        facet_wrap(~parameter, ncol = col_count, nrow = row_count) +
        geom_boxplot(size = 2) +
        geom_jitter(size = 3, width = 0.05) +
        # stat_compare_means(comparisons = comparisons, method = 't.test', size = 8, aes(label=..p.adj..), p.adjust.method = "bonferroni") +
        stat_compare_means(comparisons = comparisons, method = 't.test', size = 8) +
        theme_cowplot(font_family = 'Arial') + 
        xlab(paste0('SNP ', snpID, ' genotype')) +
        ylab(paste0('Parameter coefficient\nlog10(probability of deletion)')) +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        theme(legend.position = 'none', text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20)) 
    
    ggsave(file_name, plot = plot, width = 6*col_count, height = 9*row_count, units = 'in', dpi = 750, device = cairo_pdf)
}
