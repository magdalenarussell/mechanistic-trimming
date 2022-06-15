source(paste0('plotting_scripts/model_group_functions/', MODEL_GROUP, '.R'))
source(paste0(PROJECT_PATH, '/plotting_scripts/plot_paths.R'))

set_color_palette <- function(model_type_list){
    require(RColorBrewer)
    model_types = model_type_list[!(model_type_list %in% c('null', '2x4motif'))]
    colors = c(brewer.pal(7, 'Dark2'), brewer.pal(8, 'Set1'), brewer.pal(7, 'Set2'))
    colors = colors[!(colors %in% c("#E41A1C", "#FFFF33", "#FFD92F"))]
    names(colors) = model_types
    temp = c("#E41A1C", "#666666")
    names(temp) = c('2x4motif', 'null')
    colors = c(colors, temp)
    return(colors)
}

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
        position = str_remove(position, '_std')
        if (position %like% '_5end_'){
            val = -1 * as.numeric(substring(position, nchar(position), nchar(position)))
        } else if (position %like% '_3end_'){
            val = as.numeric(substring(position, nchar(position), nchar(position)))
        } else if (position %like% '_trim_'){
            val = 0
        }
        values = c(values, val)
    }
    together = data.table(positions = positions, values = values)
    return(together)
}

plot_base_count_coefficient_heatmap_single_group <- function(model_coef_matrix, file_name, with_values = FALSE, limits = NULL, write_plot = TRUE){
    model_coef_matrix = model_coef_matrix[(parameter %like% 'base_count') | (parameter %like% 'dinuc_count')]
    if (is.na(unique(model_coef_matrix$base)) | unique(model_coef_matrix$base) == ''){
        model_coef_matrix[, base := sapply(parameter, function(x) str_split(x, '_count_')[[1]][2])]
        model_coef_matrix[, parameter := sapply(parameter, function(x) str_split(x, '_count_')[[1]][1])]
        model_coef_matrix[, parameter := paste0(parameter, '_count')]
    }

    model_coef_matrix[, parameter := str_replace(parameter, '_base', '')]
    model_coef_matrix[, parameter := str_replace(parameter, '_dinuc', '')]

    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$coefficient/log(10)
    
    # order variables
    left_vars = unique(model_coef_matrix[parameter %like% 'left']$parameter)
    right_vars = unique(model_coef_matrix[parameter %like% 'right']$parameter)
    unique_bases = unique(model_coef_matrix$base)
    if (length(left_vars) == 0){
        fake_row = data.table(coefficient = rep(0, length(unique_bases)), parameter = rep("left_count", length(unique_bases)), base = unique_bases, model_group = rep(MODEL_GROUP,length(unique_bases)), log_10_pdel = rep(0,length(unique_bases)))
        extended_data = rbind(model_coef_matrix, fake_row, fill = TRUE)
    } else if (length(right_vars) == 0){
        fake_row = data.table(coefficient = rep(0, length(unique_bases)), parameter = rep("right_count", length(unique_bases)), base = unique_bases, model_group = rep(MODEL_GROUP,length(unique_bases)), log_10_pdel = rep(0,length(unique_bases)))
        extended_data = rbind(model_coef_matrix, fake_row, fill = TRUE)
    } else {
        extended_data = model_coef_matrix
    }
 
    vars = c(unique(extended_data[parameter %like% 'left']$parameter), unique(extended_data[parameter %like% 'right']$parameter))
    extended_data$parameter = factor(extended_data$parameter, levels = vars)

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    plot = ggplot(extended_data, aes(x=parameter, y=base, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Base count') +
        ylab('Base type') +
        geom_vline(xintercept = 1 + 0.5, size = 3, color = 'black') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits) +
        annotate("text", x = 0.55, y = 0.45, label = "5\'", size = 6) +  
        annotate("text", x = 2.5 , y = 0.45, label = "3\'", size = 6) 
    
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = model_coef_matrix, aes(x = parameter, y = base, label = round(log_10_pdel, 2)))
    }

    if (isTRUE(write_plot)){
        ggsave(file_name, plot = plot, width = 8, height = 4, units = 'in', dpi = 750, device = cairo_pdf)
    } else {
        return(plot)
    }
}


plot_melting_coefficient_heatmap_single_group <- function(model_coef_matrix, file_name, with_values = FALSE, limits = NULL, write_plot = TRUE){
    model_coef_matrix = model_coef_matrix[parameter %like% 'terminal_melting']
    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$coefficient/log(10)
    # order variables
    vars = c(unique(model_coef_matrix[parameter %like% 'left_terminal_melting']$parameter), unique(model_coef_matrix[parameter %like% 'right_terminal_melting']$parameter))
    model_coef_matrix$parameter = factor(model_coef_matrix$parameter, levels = vars)

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
        ggsave(file_name, plot = plot, width = 8, height = 3, units = 'in', dpi = 750, device = cairo_pdf)
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

plot_shape_coefficient_heatmap_single_group <- function(model_coef_matrix, file_name, with_values = FALSE, limits = NULL, write_plot = TRUE){
    model_coef_matrix = model_coef_matrix[(parameter %like% 'EP') | (parameter %like% 'MGW') | (parameter %like% 'HelT') | (parameter %like% 'Roll') | (parameter %like% 'ProT')]

    position_values = map_positions_to_values(unique(model_coef_matrix$parameter))
    together = merge(model_coef_matrix, position_values, by.x = 'parameter', by.y = 'positions')

    # convert to log_10
    together$log_10_pdel = together$coefficient/log(10)
    together[, shape := sapply(together$parameter, function(x) str_split(x, '_')[[1]][1])]
    # order variables
    together$values = factor(together$values)

    if (is.null(limits)){
        max_val = max(abs(together$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    motif_length = RIGHT_NUC_MOTIF_COUNT + LEFT_NUC_MOTIF_COUNT

    # base data
    base_data = together[shape %in% c('EP', 'MGW', 'ProT')]
    bond_data = together[shape %in% c('Roll', 'HelT')]

    base_plot = ggplot(base_data, aes(x=values, y=base, fill=log_10_pdel)) +
        facet_wrap(~ shape, ncol = 1, strip.position='left') +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Base position') +
        ylab('') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank()) +
        geom_vline(xintercept = LEFT_NUC_MOTIF_COUNT + 0.5, size = 3, color = 'black') +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits) +
        annotate("text", x = 0.45, y = 0.55, label = "5\'", size = 6) +  
        annotate("text", x = motif_length + 0.55, y = 0.55, label = "3\'", size = 6) 
 
   bond_plot = ggplot(bond_data, aes(x=values, y=base, fill=log_10_pdel)) +
        facet_wrap(~ shape, ncol = 1, strip.position='left') +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Bond position') +
        ylab('') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank()) +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits) +
        annotate("text", x = 0.45, y = 0.55, label = "5\'", size = 6) +  
        annotate("text", x = motif_length-1 + 0.55, y = 0.55, label = "3\'", size = 6) 
    
    if (with_values == TRUE){
        bond_plot = bond_plot +
            geom_text(data = bond_data, aes(x = values, y = base, label = round(log_10_pdel, 2)))
        base_plot = base_plot +
            geom_text(data = base_data, aes(x = values, y = base, label = round(log_10_pdel, 2)))
    }

    # align plots
    legend = get_legend(base_plot + theme(legend.box.margin = margin(0, 0, 0, 12)))
    plots_tog = plot_grid(base_plot + theme(legend.position = 'none'), bond_plot + theme(legend.position = 'none'), ncol = 1, rel_heights = c(3, 2))
    all_tog = plot_grid(plots_tog, legend, rel_widths = c(1, 0.9))

    if (isTRUE(write_plot)){
        ggsave(file_name, plot = all_tog, width = motif_length + 7, height = 6, units = 'in', dpi = 750, device = cairo_pdf)
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
    if (model_type %like% 'two_side_terminal') {
        model_type_name = paste0(model_type, '_', terminal_melting_5_end_length_filter)
    } else {
        model_type_name = model_type
    }
    file_name = paste0(path, '/model_evaluation_heatmap_model_type_filter_', model_type_name, '.pdf')
    eval_data = process_model_evaluation_file(eval_data, model_types_neat = model_type, terminal_melting_5_end_length_filter = terminal_melting_5_end_length_filter)
    eval_data = eval_data[motif_type == MOTIF_TYPE]

    fill_lab = make_loss_type_names_neat(type)

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

plot_model_evaluation_hairpin_nick_paracoord <- function(all_eval_data, model_type_list, type_of_loss, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, custom_name = NULL, loss_bound = NULL, color_palette = NULL, plot_size = NULL) {
    # process evaluation file and combine with Murugan 2x4 motif evaluation losses
    eval_data_murugan = process_model_evaluation_file(all_eval_data, 'motif', 2, 4, NA)
    eval_data_murugan$model_type = mapvalues(eval_data_murugan$model_type, from = 'motif', to = '2x4motif')
    eval_data = process_model_evaluation_file(all_eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter)
    eval_data = rbind(eval_data, eval_data_murugan)
    nice_names = make_model_names_neat(unique(eval_data$model_type)) 
    eval_data$nice_model_type = mapvalues(eval_data$model_type, from = unique(eval_data$model_type), to = nice_names)

    nice_hairpin_names = make_hairpin_names_neat(unique(eval_data$motif_type)) 
    eval_data$nice_nick_type = mapvalues(eval_data$motif_type, from = unique(eval_data$motif_type), to = nice_hairpin_names)
    eval_data$nice_nick_type = str_replace_all(eval_data$nice_nick_type, 'opening', '\n')
    
    eval_data$nice_nick_type = factor(eval_data$nice_nick_type, levels = c('-2 hairpin \n position', '-1 hairpin \n position', 'blunt hairpin \n position', '+1 hairpin \n position', '+2 hairpin \n position', '+3 hairpin \n position'))

    # reformat loss names to be nice
    eval_data = eval_data[loss_type == type_of_loss]
    eval_data$nice_loss_type = make_loss_type_names_neat(type_of_loss) 

    # create label dataset
    label_data = eval_data[nice_nick_type == '+3 hairpin \n position'] 

    
    # create plot
    require(ggrepel)
    plot = ggplot(eval_data) +
        geom_point(aes(y = loss, x = nice_nick_type, color = nice_model_type), size = 14)+
        geom_line(aes(y = loss, x = nice_nick_type, group = nice_model_type, color = nice_model_type), size = 9, alpha = 0.8)+
        geom_text_repel(data = label_data, aes(y = loss, x = nice_nick_type, label = nice_model_type, color = nice_model_type), nudge_x = 0.2, fontface = "bold", size = 12, direction = 'y', hjust = 0, point.padding = 1, max.overlaps = Inf) +
        theme_cowplot(font_family = 'Arial') + 
        xlab(' ') +
        ylab('Log loss\n') +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        theme(legend.position = 'none', text = element_text(size = 36), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 35), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))  +
        scale_x_discrete(expand = expansion(add = c(0.2, 4)))

    if (!is.null(loss_bound)){
        plot = plot +
            ylim(loss_bound)
    }

    if (!is.null(color_palette)){
        plot = plot + scale_color_manual(values = color_palette)
    }

    path = get_model_eval_file_path(type_of_loss)
    file_name = paste0(path, '/neat_loss_paracoord_', left_motif_size_filter, '_', right_motif_size_filter)
    if (!is.null(custom_name)){
        file_name = paste0(file_name, '_', custom_name, '.pdf')
    } else {
        file_name = paste0(file_name, '_hairpin_nick_position_compare.pdf')
    }

    if (!is.null(plot_size)){
        width = plot_size[1]
        height = plot_size[2]
    } else {
        width = 32
        height = 25
    }

    ggsave(file_name, plot = plot, width = width, height = height, units = 'in', dpi = 750, device = cairo_pdf)
}


plot_model_evaluation_loss_paracoord <- function(all_eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, custom_name = NULL, loss_bound = NULL, color_palette = NULL, same_motif_type = TRUE, plot_size = NULL) {
    # process evaluation file and combine with Murugan 2x4 motif evaluation losses
    eval_data_murugan = process_model_evaluation_file(all_eval_data, 'motif', 2, 4, NA)
    eval_data_murugan$model_type = mapvalues(eval_data_murugan$model_type, from = 'motif', to = '2x4motif')
    eval_data = process_model_evaluation_file(all_eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter)
    eval_data = rbind(eval_data, eval_data_murugan)
    if (isTRUE(same_motif_type)){
        eval_data = eval_data[motif_type == MOTIF_TYPE]
        # reformat model names to be nice
        nice_names = make_model_names_neat(unique(eval_data$model_type)) 
        eval_data$nice_model_type = mapvalues(eval_data$model_type, from = unique(eval_data$model_type), to = nice_names)
    } else {
        stopifnot(model_type_list == MODEL_TYPE)
        eval_data = eval_data[model_type == MODEL_TYPE]
        nice_hairpin_names = make_hairpin_names_neat(unique(eval_data$motif_type)) 
        eval_data$nice_model_type = mapvalues(eval_data$motif_type, from = unique(eval_data$motif_type), to = nice_hairpin_names)
    }

    # reformat loss names to be nice
    nice_loss_names = c()
    for (type in unique(eval_data$loss_type)){
        nice_name = make_loss_type_names_neat(type)
        nice_loss_names = c(nice_loss_names,  nice_name)
    }
    eval_data$nice_loss_type = mapvalues(eval_data$loss_type, from = unique(eval_data$loss_type), to = nice_loss_names)
    
    # create label dataset
    ordered_losses = order_losses(nice_loss_names)
    eval_data$nice_loss_type = factor(eval_data$nice_loss_type, levels = ordered_losses)
    last_loss = ordered_losses[length(ordered_losses)]
    label_data = eval_data[nice_loss_type == last_loss] 

    # create plot
    require(ggrepel)
    plot = ggplot(eval_data) +
        geom_point(aes(y = loss, x = nice_loss_type, color = nice_model_type), size = 14)+
        geom_line(aes(y = loss, x = nice_loss_type, group = nice_model_type, color = nice_model_type), size = 9, alpha = 0.8)+
        geom_text_repel(data = label_data, aes(y = loss, x = nice_loss_type, label = nice_model_type, color = nice_model_type), nudge_x = 0.2, fontface = "bold", size = 12, direction = 'y', hjust = 0, point.padding = 1, max.overlaps = Inf) +
        theme_cowplot(font_family = 'Arial') + 
        xlab(' ') +
        ylab('Log loss\n') +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        theme(legend.position = 'none', text = element_text(size = 36), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 35), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))  +
        scale_x_discrete(expand = expansion(add = c(0.2, 4)))

    if (!is.null(loss_bound)){
        plot = plot +
            ylim(loss_bound)
    }

    if (!is.null(color_palette)){
        plot = plot + scale_color_manual(values = color_palette)
    }

    path = get_model_eval_file_path('compare')
    file_name = paste0(path, '/neat_loss_paracoord_', left_motif_size_filter, '_', right_motif_size_filter)
    if (!is.null(custom_name)){
        file_name = paste0(file_name, '_', custom_name, '.pdf')
    } else {
        file_name = paste0(file_name, '_motif.pdf')
    }

    if (!is.null(plot_size)){
        width = plot_size[1]
        height = plot_size[2]
    } else {
        width = 48 
        height = 25
    }

    ggsave(file_name, plot = plot, width = width, height = height, units = 'in', dpi = 750, device = cairo_pdf)
}


plot_model_evaluation_scatter_coef_count <- function(eval_data, type, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, color_palette = NULL) {
    # process evaluation file
    eval_data$loss_type = type
    setnames(eval_data, type, 'loss')
    eval_data_murugan = process_model_evaluation_file(eval_data, 'motif', 2, 4, NA)
    eval_data_murugan$model_type = mapvalues(eval_data_murugan$model_type, from = 'motif', to = '2x4motif')
    eval_data = process_model_evaluation_file(eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter)

    eval_data = rbind(eval_data, eval_data_murugan)
    eval_data = eval_data[motif_type == MOTIF_TYPE]

    ylab = make_loss_type_names_neat(type)
    
    nice_names = make_model_names_neat(unique(eval_data$model_type)) 
    eval_data$nice_model_type = mapvalues(eval_data$model_type, from = unique(eval_data$model_type), to = nice_names)

    label_cols = c('loss', 'model_parameter_count', 'nice_model_type', 'held_out_clusters')
    label_data = unique(eval_data[, ..label_cols])

    require(ggrepel)
    plot = ggplot(eval_data) +
        geom_point(aes(y = loss, x = model_parameter_count, color = nice_model_type), size = 5)+
        theme_cowplot(font_family = 'Arial') + 
        xlab('Total number of terms') +
        ylab(ylab) +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        ylim(min(eval_data$loss)-0.05, max(eval_data$loss)+0.05) +
        coord_cartesian(clip = 'off') +
        geom_text_repel(data = label_data, aes(y = loss, x = model_parameter_count, label = nice_model_type, color = nice_model_type), nudge_x = 0.75, hjust = 0, point.padding = 1, max.overlaps = Inf, fontface = "bold", size = 6, xlim = c(0, 90), ylim = c(min(eval_data$loss)-0.05, max(eval_data$loss)+0.05), direction = 'y') +
        theme(legend.position = 'none', text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 18), plot.margin = unit(c(0.5,12,0.5,0.5), "cm")) 

    if (type == 'v_gene_family_loss') {
        plot = plot + facet_grid(cols = vars(held_out_clusters))
        width_dim = 12 * length(unique(eval_data$held_out_clusters)) 
    } else {
        width_dim = 16
    }

    path = get_model_eval_file_path(type)
    file_name = paste0(path, '/neat_', type, '_term_count_scatter_', left_motif_size_filter, '_', right_motif_size_filter, '_motif.pdf')

    if (!is.null(color_palette)){
        plot = plot + scale_color_manual(values = color_palette)
    }

    ggsave(file_name, plot = plot, width = width_dim, height = 8.5, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_model_evaluation_compare <- function(eval_data1, eval_data2, type1, type2, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, label = FALSE) {
    # merge two types
    cols = colnames(eval_data1)
    cols = cols[!(cols %in% c('held_out_gene_fraction', 'sample_repetitions'))]
    eval_data = merge(eval_data1, eval_data2, by = cols[!(cols == type1)])

    # process evaluation file
    eval_data = process_model_evaluation_file(eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter)
    eval_data = eval_data[motif_type == MOTIF_TYPE]

    xlab = make_loss_type_names_neat(type1)
    ylab = make_loss_type_names_neat(type2)
    
    nice_names = make_model_names_neat(unique(eval_data$model_type)) 
    eval_data$nice_model_type = mapvalues(eval_data$model_type, from = unique(eval_data$model_type), to = nice_names)
 
    plot = ggplot(eval_data) +
        geom_abline(intercept = 0, slope = 1, size = 3, color = 'gray60', linetype = 'dashed') +
        geom_point(aes(y = get(type2), x = get(type1), color = nice_model_type), size = 5)+
        theme_cowplot(font_family = 'Arial') + 
        xlab(xlab) +
        ylab(ylab) +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 

    path = get_model_eval_file_path('compare')
    file_name = paste0(path, '/neat_', type1, '_', type2, '_eval_compare_', left_motif_size_filter, '_', right_motif_size_filter, '_motif')

    if (isTRUE(label)){
        require(ggrepel)
        plot = plot +
            coord_cartesian(clip = 'off') +
            geom_text_repel(aes(y = get(type2), x = get(type1), label = nice_model_type, color = nice_model_type), nudge_x = 0.01, hjust = 0, point.padding = 1, max.overlaps = Inf, fontface = "bold", size = 6, xlim = c(0, 5), direction = 'y') +
            theme(legend.position = 'none', text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 18), plot.margin = unit(c(0.5,15,0.5,0.5), "cm")) 

        file_name = paste0(file_name, '.pdf')
        ggsave(file_name, plot = plot, width = 15, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
    } else {
        plot = plot +
            theme(text = element_text(size = 25), axis.line = element_blank(), axis.ticks = element_blank()) 
        file_name = paste0(file_name, '_no_label.pdf')
        ggsave(file_name, plot = plot, width = 18, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
    }
}

plot_model_evaluation_compare_bounds <- function(eval_data, upper_bound1, upper_bound2, type, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, label = FALSE, limits = NULL) {
    # process evaluation file
    eval_data1 = process_model_evaluation_file(eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, upper_trim_bound = upper_bound1)
    eval_data2 = process_model_evaluation_file(eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, upper_trim_bound = upper_bound2)

    cols = colnames(eval_data1)
    together = merge(eval_data1, eval_data2, by = cols[!(cols %in% c('upper_bound', 'model_parameter_count', type))])
    together = together[motif_type == MOTIF_TYPE]

    lab = make_loss_type_names_neat(type)
 
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

plot_coefficient_by_snp <- function(coef_snp_data, snpID, parameter_group = NULL, negative_control = FALSE, negative_control_type = NULL){
    require(ggpubr)
    file_path = get_individual_comparison_file_path(negative_control, negative_control_type)
    if (!is.null(parameter_group)){
        coef_snp_data = coef_snp_data[parameter %like% parameter_group]
        file_name = paste0(file_path, '/', parameter_group, '_parameters_by_snp', snpID, '.pdf')
        if (unique(coef_snp_data$parameter) == 'as.factor(trim_length)'){
            coef_snp_data$parameter = factor(coef_snp_data$parameter, levels = paste0('trim_length_', seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)))
        } else if (parameter_group == 'motif'){
            motif_order = paste0('motif_', c(rep(5, LEFT_NUC_MOTIF_COUNT), rep(3, RIGHT_NUC_MOTIF_COUNT)), 'end_pos', c(seq(LEFT_NUC_MOTIF_COUNT, 1), seq(1, RIGHT_NUC_MOTIF_COUNT)), ', base ', rep(c('A', 'C', 'G', 'T'), each = LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT))
            coef_snp_data$parameter = factor(coef_snp_data$parameter, levels = motif_order)
        }
    } else {
        file_name = paste0(file_path, '/parameters_by_snp', snpID, '.pdf')
    }
    subset = coef_snp_data[snp == snpID & !is.na(genotype)]
    subset$genotype = as.character(subset$genotype)
    subset$log_10_coef = subset$coefficient/log(10)
    subset = subset[!is.na(parameter)]
    comparisons = list(c("0", "1"), c("1", "2"), c("0", "2")) 
    if (parameter_group == 'motif'){
        row_count = 4
        col_count = LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT
    } else {
        row_count = ceiling(length(unique(subset$parameter))/6)
        col_count = min(6, length(unique(subset$parameter)))
    }
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
