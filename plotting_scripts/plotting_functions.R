source(paste0('plotting_scripts/model_group_functions/', MODEL_GROUP, '.R'))
source(paste0(PROJECT_PATH, '/scripts/data_compilation_functions.R'))

compile_trims <- function(directory){
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    together = foreach(file = files, .combine = 'rbind') %dopar% {
        temp = fread(file)
        subject_id = extract_subject_ID(file)
        temp$subject = subject_id
        cols = c('subject', 'productive', GENE_NAME, TRIM_TYPE)
        temp[, .N, by = cols]
    }
    stopImplicitCluster()
    
    new_cols = c('subject', 'productive', TRIM_TYPE)
    concat = together[, sum(N), by = new_cols]
    setnames(concat, 'V1', 'count')
    return(concat)
}

get_gene_sequence <- function(gene_name, gene_seq_length, pnuc_count = 2){
    whole_nucseq = get_oriented_whole_nucseqs()
    temp_data = whole_nucseq[substring(gene, 4, 4) == substring(GENE_NAME, 4, 4)]
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

get_predicted_dist_figure_file_path <- function(){
    if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'predicted_trimming_distributions')
    dir.create(path, recursive = TRUE)
    return(path)
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

get_coef_heatmap_file_path <- function(){
    if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
    } else {
        model = MODEL_TYPE
    }

    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'model_coefficient_heatmaps')
    dir.create(path, recursive = TRUE)
    return(path)
}


plot_model_coefficient_heatmap_single_group <- function(model_coef_matrix, file_name, with_values = FALSE, limits = NULL, write_plot = TRUE){
    data = model_coef_matrix %>%
        pivot_longer(!c(base, model_group), names_to = 'position', values_to = 'coef')

    position_values = map_positions_to_values(unique(data$position))
    together = merge(data, position_values, by.x = 'position', by.y = 'positions')

    # convert to log_10
    together$log_10_pdel = together$coef/log(10)
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
        # scale_fill_viridis_c(name = 'log10(probability of deletion)', limits = limits)
    
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = together, aes(x = values, y = base, label = round(log_10_pdel, 2)))
    }

    if (isTRUE(write_plot)){
        ggsave(file_name, plot = plot, width = 9, height = 4, units = 'in', dpi = 750, device = cairo_pdf)
    } else {
        return(plot)
    }
}

get_model_coef_heatmap_compare_file_name <- function(model_type1, model_type2, left_side_terminal_melt_length1 = NA, left_side_terminal_melt_length2 = NA){
    if (grepl('_side_terminal_melting', model_type1, fixed = TRUE)){
        stopifnot(!(is.na(left_side_terminal_melt_length1)))
        model1 = paste0(model_type1, '_', left_side_terminal_melt_length1, '_length_melting_left')
    } else {
        model1 = model_type1
    }

    if (grepl('_side_terminal_melting', model_type2, fixed = TRUE)){
        stopifnot(!(is.na(left_side_terminal_melt_length2)))
        model2 = paste0(model_type2, '_', left_side_terminal_melt_length2, '_length_melting_left')
    } else {
        model2 = model_type2
    }

    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, MODEL_GROUP, GENE_WEIGHT_TYPE, paste0('compare_', model1, '-', model2) , paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
    dir.create(path, recursive = TRUE)
    
    name = file.path(path, 'heatmap_compare.pdf')
    return(name)
}

plot_model_coefficient_heatmap_compare <- function(model_coef_matrix1, model_coef_matrix2, file_name, with_values = FALSE, limits = NULL, write_plot = TRUE){
    pos = get_positions()
    #get difference between two matrices
    data = model_coef_matrix1[model_coef_matrix2, (pos) := Map(`/`, mget(pos), mget(paste0("i.", pos))), on = .(model_group, base)]
    data = data %>%
        pivot_longer(!c(base, model_group), names_to = 'position', values_to = 'coef')

    position_values = map_positions_to_values(unique(data$position))
    together = merge(data, position_values, by.x = 'position', by.y = 'positions')

    # convert to log_10
    # order variables
    together$base = factor(together$base, levels = c('T', 'G', 'C', 'A'))
    together$values = factor(together$values)

    if (is.null(limits)){
        max_val = max(abs(together$coef))
        limits = c(-max_val, max_val)
    }
    
    plot = ggplot(together, aes(x=values, y=base, fill=coef)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Position') +
        ylab ('Base') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        geom_vline(xintercept = RIGHT_NUC_MOTIF_COUNT + 0.5, size = 3, color = 'black') +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits) +
        annotate("text", x = 0.35, y = 0.32, label = "5\'", size = 8) +  
        annotate("text", x = 8.65, y = 0.32, label = "3\'", size = 8) +  
        coord_cartesian(ylim = c(1, 4), clip = "off")
        # scale_fill_viridis_c(name = 'log10(probability of deletion)', limits = limits)
    
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = together, aes(x = values, y = base, label = round(coef, 2)))
    }

    if (isTRUE(write_plot)){
        ggsave(file_name, plot = plot, width = 9, height = 4, units = 'in', dpi = 750, device = cairo_pdf)
    } else {
        return(plot)
    }
}



get_residual_figure_file_path <- function(){
    if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
    } else {
        model = MODEL_TYPE
    }

    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'predicted_trimming_residuals')

    dir.create(path, recursive = TRUE)
    return(path)
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

get_coef_variations_file_path <- function(){
    if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
    } else {
        model = MODEL_TYPE
    }

    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'model_coefficient_variations')

    dir.create(path, recursive = TRUE)
    return(path)
}

get_coef_variations_file_name <- function(){
    filename = paste0('coefficient_variations.pdf')
    return(filename)
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
    file = get_coef_variations_file_name()

    file_name = file.path(path, file)
    ggsave(file_name, plot = plot, width = 15, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
}

get_all_residual_figure_file_path <- function(){
    if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
    } else {
        model = MODEL_TYPE
    }

    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'ALL_predicted_trimming_residuals')
    dir.create(path, recursive = TRUE)
    return(path)
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

get_base_composition_file_path <- function(){
    if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
    } else {
        model = MODEL_TYPE
    }

    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'model_base_composition')
    dir.create(path, recursive = TRUE)
    return(path)
}

get_gene_composition_data <- function(motif_data, weighting = 'uniform'){
    positions = get_positions()
    if (weighting == 'uniform'){
        cols = c('gene', 'trim_length', positions)
        temp = unique(motif_data[, ..cols]) 
        temp$count = 1
    } else if (weighting == 'raw_count'){
        cols = c('gene', 'trim_length', positions)
        temp = motif_data[, sum(count), by = cols]
        setnames(temp, 'V1', 'count')
    } else if (weighting == 'p_gene_marginal'){
        cols = c('gene', 'trim_length', 'p_gene', positions)
        temp = unique(motif_data[, ..cols]) 
        temp$count = temp$p_gene
    } else if (weighting == 'weighted_observation'){
        temp = motif_data
        temp$count = temp$weighted_observation
    }

    composition = temp[, sum(count), by = positions]
    setnames(composition, 'V1', 'weighted_obs_sum')

    composition_pivot = composition %>%
        pivot_longer(!c(weighted_obs_sum), names_to = 'position', values_to = 'base') %>%
        as.data.table()

    comp_sum = composition_pivot[, sum(weighted_obs_sum), by = .(position, base)]
    setnames(comp_sum, 'V1', 'weighted_obs_sum')

    return(comp_sum)
}

plot_gene_composition <- function(motif_data, weighting = 'uniform'){
    stopifnot(weighting %in% c('uniform', 'p_gene_marginal', 'weighted_observation', 'raw_count'))
    path = get_base_composition_file_path()
    file_name = paste0(path, '/base_composition_', weighting, '_weighting.pdf')

    comp_sum = get_gene_composition_data(motif_data, weighting) 

    position_values = map_positions_to_values(unique(comp_sum$position))
    together = merge(comp_sum, position_values, by.x = 'position', by.y = 'positions')
    together$values = factor(together$values)

    plot = ggplot(together, aes(values, weighted_obs_sum, fill = base)) +
        geom_bar(stat = 'identity', position = 'fill') +
        theme_cowplot(font_family = 'Arial') +
        background_grid(major = 'y') +
        xlab('Position') +
        ylab('Proportion') +
        scale_fill_brewer(palette = 'Set2') +
        theme(text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.ticks = element_line(color = 'gray60', size = 1.5)) 

    ggsave(file_name, plot = plot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
    return(plot) 
}

get_resid_compare_file_path <- function(){
    if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
    } else {
        model = MODEL_TYPE
    }

    path = file.path(PROJECT_PATH, 'plots',ANNOTATION_TYPE, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'residual_comparison')
    dir.create(path, recursive = TRUE)
    return(path)
}


plot_residual_scatter <- function(residual_mse_df, feature_df, annotate = TRUE){
    together = merge(residual_mse_df, feature_df)
    path = get_resid_compare_file_path()
    file_name = paste0(path, '/model_residual_', RESIDUAL_COMPARE_FEATURE, '_scatter.pdf')

    plot = ggplot(together, aes(y = rmse, x = feature)) +
        geom_point(size = 3, alpha = 0.6) +
        stat_cor(aes(label=..rr.label..), label.x.npc = 'left', label.y.npc = 'top') +
        theme_cowplot(font_family = 'Arial') +
        background_grid(major = 'xy') +
        ylab('Root mean square error') +
        xlab(RESIDUAL_COMPARE_FEATURE) +
        theme(text = element_text(size = 20), axis.ticks = element_line(color = 'gray60', size = 1.5)) 

    if (isTRUE(annotate)){
        plot = plot + 
            geom_text(aes(label=ifelse(rmse > quantile(together$rmse, 0.5) & feature > quantile(together$feature, 0.5), gene, '')), hjust = 0, nudge_x = 0.002)
    }

    ggsave(file_name, plot = plot, width = 10, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

    return(plot)
}


plot_positional_residual_scatter <- function(residual_avg_df, features_df, annotate = TRUE){
    path = get_resid_compare_file_path()
    file_name = paste0(path, '/model_residual_positional_slope_', RESIDUAL_COMPARE_FEATURE, '_scatter.pdf')

    model = lm(avg_resid ~ trim_length*as.factor(gene), residual_avg_df)
    residual_avg_df$resid_slope = coef(model)['trim_length'] + coef(model)[paste0('trim_length:as.factor(gene)', residual_avg_df$gene)] 
    residual_avg_df[gene == 'TRBV1', resid_slope := coef(model)['trim_length']]
    
    simple_resid_slopes = unique(residual_avg_df[, c('gene', 'resid_slope')])
    together = merge(simple_resid_slopes, features_df)

    plot = ggplot(together, aes(y = resid_slope, x = feature))+
        geom_point(size = 3, alpha = 0.6) +
        stat_cor(aes(label=..rr.label..), label.x.npc = 'left', label.y.npc = 'top') +
        theme_cowplot(font_family = 'Arial') +
        background_grid(major = 'xy') +
        ylab('Slope of average positional residuals') +
        xlab(RESIDUAL_COMPARE_FEATURE) +
        theme(text = element_text(size = 20), axis.ticks = element_line(color = 'gray60', size = 1.5)) 
    if (isTRUE(annotate)){
        plot = plot + 
            geom_text(aes(label=ifelse(resid_slope < quantile(together$resid_slope, 0.5) & feature > quantile(together$feature, 0.5), gene, '')), hjust = 0, nudge_x = 0.002)
    }

    ggsave(file_name, plot = plot, width = 10, height = 10, units = 'in', dpi = 750, device = cairo_pdf)

    return(plot)
}



get_model_eval_file_path <- function(type){
    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, MODEL_GROUP, paste0('model_evaluation_', type), paste0(TRIM_TYPE, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
    dir.create(path, recursive = TRUE)
    return(path)
}

plot_model_evaluation_heatmap <- function(type, with_values = FALSE, model_type_filter = NA){
    file_path = get_model_evaluation_file_name(type) 
    path = get_model_eval_file_path(type)
    file_name = paste0(path, '/model_evaluation_heatmap_model_type_filter_', model_type_filter, '.pdf')
    eval_data = fread(file_path)
    eval_data = process_model_evaluation_file(eval_data, model_types_neat = model_type_filter)
    eval_data = eval_data[motif_type == MOTIF_TYPE]

    if (type == 'log_loss'){
        fill_lab = 'Conditional log loss'
    } else if (type == 'per_gene'){
        fill_lab = 'Average root mean\nsquared error\nacross genes'
    } else if (type == 'per_gene_per_trim'){
        fill_lab = 'Average root mean\nsquared error\nacross genes, trims'
    }

    setnames(eval_data, type, 'loss')

    plot = ggplot(eval_data) +
        facet_wrap(~model_type) +
        geom_tile(aes(x = motif_length_5_end, y = motif_length_3_end, fill = loss)) +
        theme_cowplot(font_family = 'Arial') + 
        xlab('5\' motif length') +
        ylab('3\' motif length') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_viridis_c(name = fill_lab)
    
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = eval, aes(x = motif_length_5_end, y = motif_length_3_end, label = round(loss, 0)))
    }

    ggsave(file_name, plot = plot, width = 10, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_model_evaluation_scatter <- function(type, model_type_list = NA){
    file_path = get_model_evaluation_file_name(type)
    path = get_model_eval_file_path(type)
    if (length(model_type_list) == 1){
        file_name = paste0(path, '/model_evaluation_scatter_model_type_filter_', model_type_list, '.pdf')
    } else {
        file_name = paste0(path, '/model_evaluation_scatter_model_type_filter_NA.pdf')
    }
    eval_data = fread(file_path)
    eval_data = process_model_evaluation_file(eval_data, model_types_neat = model_type_list)
    eval_data = eval_data[motif_type == MOTIF_TYPE]

    eval_data[, length := motif_length_5_end + motif_length_3_end]

    setnames(eval_data, type, 'loss')
    eval_data = unique(eval_data[, c('length', 'loss', 'model_type', 'motif_length_5_end')])

    if (type == 'log_loss'){
        ylab = 'Conditional log loss'
    } else if (type == 'per_gene'){
        ylab = 'Average root mean squared error across genes'
    } else if (type == 'per_gene_per_trim'){
        ylab = 'Average root mean squared error across genes, trims'
    }

    if (length(unique(eval_data$model_type)) > 3){
        plot = ggplot(eval_data) +
            geom_point(aes(x = length, y = loss, color = model_type), size = 6, alpha = 0.7) 
    } else {
        plot = ggplot(eval_data) +
            geom_point(aes(x = length, y = loss, shape = model_type, color = motif_length_5_end), size = 6, alpha = 0.7) +
            scale_color_viridis_c(name = '5\' motif length')
    }
    
    final_plot = plot +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Total motif length') +
        ylab(ylab) +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        guides(fill = guide_colourbar(barheight = 14)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 

    ggsave(file_name, plot = final_plot, width = 13, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_average_trims <- function(directory){
    data = compile_trims(directory)
    expanded = as.data.table(lapply(data, rep, data$count))
    expanded[productive == TRUE, productivity:= 'productive']
    expanded[productive == FALSE, productivity:= 'non-productive']

    plot = ggplot(expanded, aes(x=get(TRIM_TYPE), group = subject)) +
           facet_grid(rows = vars(productivity)) +
           stat_ecdf(geom = "step", size = 0.5, alpha = 0.5) +
           theme_cowplot(font_family = 'Arial') +
           theme(legend.position = "none", axis.text = element_text(size = 24), panel.spacing = unit(2, "lines"), strip.text = element_text(size = 22), axis.line = element_blank(), text = element_text(size = 40), axis.ticks = element_line(color = 'gray60', size = 1.5)) +
            background_grid(major = 'xy') +
            xlab(TRIM_TYPE) +
            ylab('Cumulative probability')

    filename = paste0(PROJECT_PATH, '/plots/random/', TRIM_TYPE, '_dist_by_subject.pdf')
    ggsave(filename, plot = plot, width = 12, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_motif_coefficient_distribution <<- function(model_coef_matrix){
    data = model_coef_matrix %>%
        pivot_longer(!c(base, model_group), names_to = 'position', values_to = 'coef') %>%
        as.data.table()
    
    data$log_10_pdel = data$coef/log(10)
   
    plot = ggplot(data, aes(x=log_10_pdel)) + 
        geom_histogram(aes(y=..density..), binwidth=.1, fill = 'blue', alpha = 0.6) +
        geom_density(alpha=.8, color="black", size = 3) +
        xlim(-1, 1)
   
    final_plot = plot +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Motif coefficient') +
        ylab('Probability density')
        ylab(ylab) +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 

    return(final_plot)
}

scatter_two_covariates <- function(motif_data, cov1, cov2){
    plot = ggplot(motif_data) + 
        geom_jitter(aes(x = get(cov1), y = get(cov2)), size = 4, alpha = 0.4)
        geom_histogram(aes(y=..density..), binwidth=.1, fill = 'blue', alpha = 0.6) +
        geom_density(alpha=.8, color="black", size = 3) +
        xlim(-1, 1)
   
    final_plot = plot +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Motif coefficient') +
        ylab('Probability density')
        ylab(ylab) +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 


}

plot_pwm_profiles <- function(condensed_pwm_data, model_name){
    max_data = condensed_pwm_data[, max_pwm_score := max(pwm_score), by = gene]
    max_data = max_data[max_pwm_score == pwm_score]
    median = condensed_pwm_data[, median(pwm_score), by = .(trim_length)]
    setnames(median, 'V1', 'median_pwm_score')
    
    plot = ggplot(condensed_pwm_data) + 
        geom_line(aes(x=trim_length, y = pwm_score, group = gene), size = 2, alpha = 0.6)+
        geom_line(data = median, mapping = aes(x=trim_length, y = median_pwm_score), size = 2, color = 'blue')+
        geom_point(data = max_data, mapping = aes(x=trim_length, y = pwm_score), size = 4, color = 'red', alpha = 0.6) +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Trimming length') +
        ylab(paste0('Motif PWM score')) +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 
    

}
