source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/model_group_functions/', MODEL_GROUP, '.R'))
source(paste0(MOD_PROJECT_PATH, '/plotting_scripts/plot_paths.R'))

set_color_palette <- function(model_type_list, with_params = FALSE){
    require(RColorBrewer)
    if (isTRUE(with_params)){
        model_types = model_type_list[!(model_type_list %in% c('null (0 params)', '2x4motif (18 params)'))]
    } else {
        model_types = model_type_list[!(model_type_list %in% c('null', '2x4motif'))]
    }
    colors = c(brewer.pal(7, 'Dark2'), brewer.pal(8, 'Set1'), brewer.pal(7, 'Set2'))
    colors = colors[!(colors %in% c("#E41A1C", "#FFFF33", "#FFD92F"))]
    names(colors) = model_types
    temp = c("#E41A1C", "#666666", "")
    if (isTRUE(with_params)){
        names(temp) = c('2x4motif (18 params)', 'null (0 params)')
    } else {
        names(temp) = c('2x4motif', 'null')
    }
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
    together = data.table(base = str_split(unlist(as.character(gene_sequence)), '')[[1]], position = positions)
    return(together)
}

get_motif_colors <- function(gene_seq_positions, motif, highlight_color){
    stopifnot(nchar(motif) == 3)
    gene_seq_positions$index = seq(1, nrow(gene_seq_positions))
    gene_seq_positions[index < 15, start_motif := base]
    gene_seq_positions[index < 15, start_motif := paste0(start_motif, gene_seq_positions$base[index + 1], gene_seq_positions$base[index + 2])] 
    for (mot in motif){
        gene_seq_positions[start_motif == mot, color := highlight_color]
        start_index = gene_seq_positions[start_motif == mot]$index
        gene_seq_positions[index == start_index + 1 | index == start_index + 2, color := highlight_color]
    }
    gene_seq_positions[color != highlight_color | is.na(color), color := 'black']
    return(gene_seq_positions)
}

plot_predicted_trimming_dists_single_group <- function(data, gene_name, file_name = NULL, ylim = NULL, write_plot = TRUE, color = 'blue', motif_highlight_color = 'black', motif_highlight = c('CTT', 'CGT'), seq_text = 7){
    important_cols = c('trim_length', 'predicted_prob', 'gene')
    predicted_data = data[gene == gene_name, ..important_cols]
    empirical_data = data[gene == gene_name][order(subject, trim_length)]
    total_tcrs = paste(data[gene == gene_name, sum(count)])
    title = paste0(gene_name, ', Total TCR count = ', total_tcrs)
    # get gene sequence
    gene_seq = get_gene_sequence(gene_name, max(data$trim_length))
    gene_seq_with_positions = get_plot_positions_for_gene_sequence(gene_seq)
    gene_seq_with_positions =get_motif_colors(gene_seq_with_positions, motif_highlight, 'motif')
        
    if (!is.null(ylim)){
        max_prob = ylim
    } else {
        max_prob = max(max(empirical_data$empirical_prob), max(predicted_data$predicted_prob))
    }

    plot = ggplot() +
        geom_line(data = empirical_data, aes(x = trim_length, y = empirical_prob, group = subject), size = 1, alpha = 0.5, color = "gray60") +
        geom_line(data = predicted_data, aes(x = trim_length, y = predicted_prob), size = 1.75, color = color) +
        geom_vline(xintercept = 0, color = 'black', size = 2) +
        geom_text(data = gene_seq_with_positions, y = max_prob, aes(x = position, label = base, color = color), size = seq_text) +
        geom_text(y = max_prob, aes(x = -2.1), label = '3\'- ', size = seq_text - 1) +
        geom_text(y = max_prob, aes(x = UPPER_TRIM_BOUND + 0.1), label = ' -5\'', size = seq_text - 1) +
        ggtitle(title) +
        xlab('Number of trimmed nucleotides') +
        ylab('Probability') +
        theme_cowplot(font_family = 'Arial') + 
        theme(legend.position = "none", text = element_text(size = 25), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        scale_color_manual(values = c(motif = motif_highlight_color, black = 'black'))
    
    if (isTRUE(write_plot)){
        stopifnot(!is.null(file_name))
        ggsave(file_name, plot = plot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
    }
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
    left_bases = unique(model_coef_matrix[parameter %like% 'left']$base)
    right_bases = unique(model_coef_matrix[parameter %like% 'right']$base)

    if (length(left_vars) == 0){
        fake_row = data.table(coefficient = rep(NA, length(unique_bases)), parameter = rep("left_count", length(unique_bases)), base = unique_bases, model_group = rep(MODEL_GROUP,length(unique_bases)), log_10_pdel = rep(NA,length(unique_bases)))
        extended_data = rbind(model_coef_matrix, fake_row, fill = TRUE)
    } else if (length(right_vars) == 0){
        fake_row = data.table(coefficient = rep(NA, length(unique_bases)), parameter = rep("right_count", length(unique_bases)), base = unique_bases, model_group = rep(MODEL_GROUP,length(unique_bases)), log_10_pdel = rep(NA,length(unique_bases)))
        extended_data = rbind(model_coef_matrix, fake_row, fill = TRUE)
    } else if (length(left_bases) < 2){
        missing = unique_bases[!(unique_bases == left_bases)]
        fake_row = data.table(coefficient = rep(NA, length(missing)), parameter = rep("left_count", length(missing)), base = missing, model_group = rep(MODEL_GROUP,length(missing)), log_10_pdel = rep(NA,length(missing)))
        extended_data = rbind(model_coef_matrix, fake_row, fill = TRUE)
    } else if (length(right_bases) < 2){
        missing = unique_bases[!(unique_bases == right_bases)]
        fake_row = data.table(coefficient = rep(NA, length(missing)), parameter = rep("right_count", length(missing)), base = missing, model_group = rep(MODEL_GROUP,length(missing)), log_10_pdel = rep(NA,length(missing)))
        extended_data = rbind(model_coef_matrix, fake_row, fill = TRUE)
    } else {
        extended_data = model_coef_matrix
    }
    require(plyr) 
    extended_data$parameter = mapvalues(extended_data$parameter, from = unique(extended_data$parameter), to = c('5\' base count', '3\' base count'))
    vars = c(unique(extended_data[parameter %like% '5']$parameter), unique(extended_data[parameter %like% '3']$parameter))
    extended_data$parameter = factor(extended_data$parameter, levels = vars)

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
    
    plot = ggplot(extended_data, aes(x=parameter, y=base, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('') +
        ylab('Base type') +
        geom_vline(xintercept = 1 + 0.5, size = 3.5, color = 'black') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits, na.value = 'gray80') 
    
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

plot_lindistance_coefficient_heatmap_single_group <- function(model_coef_matrix, file_name, with_values = FALSE, limits = NULL, write_plot = TRUE){
    model_coef_matrix = model_coef_matrix[parameter %like% 'trim_length']
    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$coefficient/log(10)
    model_coef_matrix$parameter = 'Trimming length\n(integer-valued)'

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
     
    plot = ggplot(model_coef_matrix, aes(x=parameter, y=base, fill=log_10_pdel)) +
        geom_tile() +
        geom_vline(xintercept = 0.5, size = 3.5, color = 'black') +
        theme_cowplot(font_family = 'Arial') + 
        xlab('') +
        ylab('') +
        theme(text = element_text(size = 35), axis.text.x = element_text(size = 22), axis.text.y = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
        guides(fill = guide_colourbar(barheight =12)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits) +
    
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = model_coef_matrix, aes(x = parameter, y = base, label = round(log_10_pdel, 2)), size = 10)
    }

    if (isTRUE(write_plot)){
        ggsave(file_name, plot = plot, width = 3, height = 4, units = 'in', dpi = 750, device = cairo_pdf)
    } else {
        return(plot)
    }
}

plot_distance_coefficient_heatmap_single_group <- function(model_coef_matrix, file_name, with_values = FALSE, limits = NULL, write_plot = TRUE){
    model_coef_matrix = model_coef_matrix[parameter %like% 'trim_length']
    # convert to log_10
    model_coef_matrix$log_10_pdel = model_coef_matrix$coefficient/log(10)
    # order variables
    model_coef_matrix$parameter = mapvalues(model_coef_matrix$parameter, from = paste0('trim_length_', seq(UPPER_TRIM_BOUND, LOWER_TRIM_BOUND)), to = seq(UPPER_TRIM_BOUND, LOWER_TRIM_BOUND))
    model_coef_matrix$parameter = factor(model_coef_matrix$parameter, levels = seq(UPPER_TRIM_BOUND, LOWER_TRIM_BOUND))
    total_dist = UPPER_TRIM_BOUND - LOWER_TRIM_BOUND

    if (is.null(limits)){
        max_val = max(abs(model_coef_matrix$log_10_pdel))
        limits = c(-max_val, max_val)
    }
     
    plot = ggplot(model_coef_matrix, aes(x=parameter, y=base, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('\nTrimming distance (categorical)') +
        ylab('') +
        theme(text = element_text(size = 35), axis.text.x = element_text(size = 22), axis.text.y = element_blank(), axis.line = element_blank(), axis.ticks = element_blank()) +
        guides(fill = guide_colourbar(barheight =12)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits) +
        annotate("text", x = 0.56, y = 0.45, label = "5\'", size = 6) +  
        annotate("text", x = 13.36 , y = 0.45, label = "3\'", size = 6) 
    
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
        geom_vline(xintercept = LEFT_NUC_MOTIF_COUNT + 0.5, size = 3.5, color = 'black') +
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

plot_model_evaluation_heatmap <- function(eval_data, type, with_values = FALSE, model_type, terminal_melting_5_end_length_filter, limits = NULL, write_plot = TRUE){
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
    
    if (write_plot == TRUE){
        ggsave(file_name, plot = plot, width = 11, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
    } else {
        return(plot)
    }
}

plot_model_evaluation_hairpin_nick_paracoord <- function(all_eval_data, model_type_list, type_of_loss, pre_filter = FALSE, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, custom_name = NULL, loss_bound = NULL, color_palette = NULL, plot_size = NULL, write_plot = TRUE, expand_var = 4) {
    # process evaluation file and combine with Murugan 2x4 motif evaluation losses
    if (isFALSE(pre_filter)){
        eval_data_murugan = process_model_evaluation_file(all_eval_data, 'motif', 2, 4, NA)
        eval_data_murugan$model_type = mapvalues(eval_data_murugan$model_type, from = 'motif', to = '2x4motif')
        eval_data = process_model_evaluation_file(all_eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter)
        eval_data = rbind(eval_data, eval_data_murugan)
    } else {
        eval_data = all_eval_data
    }
    
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
        scale_x_discrete(expand = expansion(add = c(0.2, expand_var)))

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
    
    if (write_plot == TRUE) {
        ggsave(file_name, plot = plot, width = width, height = height, units = 'in', dpi = 750, device = cairo_pdf)
    } else {
        return(plot)
    }
}

plot_model_evaluation_loss_paracoord <- function(all_eval_data, model_type_list, pre_filter = FALSE, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter, custom_name = NULL, loss_bound = NULL, color_palette = NULL, same_motif_type = TRUE, plot_size = NULL, write_plot = TRUE, expand_var = 4, loss_order = NULL) {
    if (isFALSE(pre_filter)){
        # process evaluation file and combine with Murugan 2x4 motif evaluation losses
        eval_data_murugan = process_model_evaluation_file(all_eval_data, 'motif', 2, 4, NA)
        eval_data_murugan$model_type = mapvalues(eval_data_murugan$model_type, from = 'motif', to = '2x4motif')
        eval_data = process_model_evaluation_file(all_eval_data, model_type_list, left_motif_size_filter, right_motif_size_filter, terminal_melting_5_end_length_filter)
        eval_data = rbind(eval_data, eval_data_murugan)
    } else {
        eval_data = all_eval_data
    }
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
    if (is.null(loss_order)){
        ordered_losses = order_losses(nice_loss_names)
    } else {
        ordered_losses = loss_order
    }
    eval_data$nice_loss_type = factor(eval_data$nice_loss_type, levels = ordered_losses)
    last_loss = ordered_losses[length(ordered_losses)]
    label_data = eval_data[nice_loss_type == last_loss] 

    # create plot
    require(ggrepel)
    plot = ggplot(eval_data) +
        geom_point(aes(y = loss, x = nice_loss_type, color = nice_model_type), size = 14)+
        geom_line(aes(y = loss, x = nice_loss_type, group = nice_model_type, color = nice_model_type), size = 9, alpha = 0.8)+
        geom_text_repel(data = label_data, aes(y = loss, x = nice_loss_type, label = nice_model_type, color = nice_model_type), nudge_x = 0.2, fontface = "bold", size = 13, direction = 'y', hjust = 0, point.padding = 1, max.overlaps = Inf, lineheight = 0.8) +
        theme_cowplot(font_family = 'Arial') + 
        xlab(' ') +
        ylab('Log loss\n') +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        theme(legend.position = 'none', text = element_text(size = 44), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 44), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))  +
        scale_x_discrete(expand = expansion(add = c(0.2, expand_var)))

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
    
    if (isTRUE(write_plot)){
        ggsave(file_name, plot = plot, width = width, height = height, units = 'in', dpi = 750, device = cairo_pdf)
    } else {
        return(plot)
    }
}
