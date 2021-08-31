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
    whole_nucseq = fread('_ignore/tcrb_processed_geneseq.tsv')
    temp_data = whole_nucseq[substring(gene_names, 4, 4) == substring(gene_name, 4, 4)]
    colnames(temp_data) = c(GENE_NAME, 'sequences')
    gene_groups = get_common_genes_from_seqs(temp_data)
    together = unique(merge(temp_data, gene_groups)[, c('gene', 'sequences')]) 
    gene = together[gene == gene_name][1]

    # get sequence
    whole_gene_seq = DNAString(gene$sequences)
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
    path = file.path(PROJECT_PATH, 'plots', 'predicted_trimming_distributions', MODEL_GROUP, GENE_WEIGHT_TYPE, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
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
            val = as.numeric(substring(position, nchar(position), nchar(position)))
        } else if (substring(position, 1, 7) == "motif_3"){
            val = -1 * as.numeric(substring(position, nchar(position), nchar(position)))
        }
        values = c(values, val)
    }
    together = data.table(positions = positions, values = values)
    return(together)
}

get_coef_heatmap_file_path <- function(){
    path = file.path(PROJECT_PATH, 'plots', 'model_coefficient_heatmaps', MODEL_GROUP, GENE_WEIGHT_TYPE, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
    dir.create(path, recursive = TRUE)

    return(path)
}


plot_model_coefficient_heatmap_single_group <- function(model_coef_matrix, file_name, with_values = FALSE, limits = NULL){
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

    plot = ggplot(together, aes(x=values, y=base, fill=log_10_pdel)) +
        geom_tile() +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Position') +
        ylab ('Base') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        geom_vline(xintercept = RIGHT_NUC_MOTIF_COUNT + 0.5, size = 3, color = 'black') +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_distiller(palette = 'PuOr', name = 'log10(probability of deletion)', limits = limits)
        # scale_fill_viridis_c(name = 'log10(probability of deletion)', limits = limits)
    
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = together, aes(x = values, y = base, label = round(log_10_pdel, 2)))
    }

    ggsave(file_name, plot = plot, width = 10, height = 5, units = 'in', dpi = 750, device = cairo_pdf)
}

get_residual_figure_file_path <- function(){
    path = file.path(PROJECT_PATH, 'plots', 'predicted_trimming_residuals', MODEL_GROUP, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
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
    path = file.path(PROJECT_PATH, 'plots', 'model_coefficient_variations', MODEL_GROUP, GENE_WEIGHT_TYPE, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
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
    path = file.path(PROJECT_PATH, 'plots', 'ALL_predicted_trimming_residuals', MODEL_GROUP, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
    dir.create(path, recursive = TRUE)
    return(path)
}


plot_all_model_residuals_plot <- function(data, file_name){
    data$residual = data$empirical_prob - data$predicted_prob
    mean_data = data[, mean(residual), by = .(gene, trim_length)]
    setnames(mean_data, 'V1', 'mean_residual')
    overall_mean_data = data[, mean(residual), by = .(trim_length)]
    setnames(overall_mean_data, 'V1', 'overall_mean_residual')
    plot = ggplot(data = mean_data, aes(x = trim_length, y = mean_residual, group = gene)) +
        geom_line(size = 2, alpha = 0.5) +
        geom_line(data = overall_mean_data, aes(x = trim_length, y = overall_mean_residual), color = 'blue', size = 2, alpha = 0.8) +
        geom_hline(yintercept = 0, color = 'gray', size = 3) +
        xlab('Number of trimmed nucleotides') +
        ylab('Observed prob - Predicted prob') +
        theme_cowplot(font_family = 'Arial') + 
        theme(legend.position = "none", text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        ylim(-0.5, 1)
    
    ggsave(file_name, plot = plot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
    return(plot) 
}

get_base_composition_file_path <- function(){
    path = file.path(PROJECT_PATH, 'plots', 'model_base_composition', MODEL_GROUP, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
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
    path = file.path(PROJECT_PATH, 'plots', 'residual_comparison', MODEL_GROUP, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
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



get_model_eval_file_path <- function(){
    path = file.path(PROJECT_PATH, 'plots', 'model_evaluation', MODEL_GROUP, paste0(TRIM_TYPE, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
    dir.create(path, recursive = TRUE)
    return(path)
}

plot_model_evaluation_heatmap <- function(file_path = get_model_evaluation_file_name(), with_values = FALSE){
    path = get_model_eval_file_path()
    file_name = paste0(path, '/model_evaluation_heatmap.pdf')
    eval_data = fread(file_path)
    plot = ggplot(eval) +
        geom_tile(aes(x = motif_length_5_end, y = motif_length_3_end, fill = log_loss)) +
        theme_cowplot(font_family = 'Arial') + 
        xlab('5\' motif length') +
        ylab('3\' motif length') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        guides(fill = guide_colourbar(barheight = 14)) +
        scale_fill_viridis_c(name = 'Conditional Log Loss')
    
    if (with_values == TRUE){
        plot = plot +
            geom_text(data = eval, aes(x = motif_length_5_end, y = motif_length_3_end, label = round(log_loss, 0)))
    }

    ggsave(file_name, plot = plot, width = 10, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
}

plot_model_evaluation_scatter <- function(file_path = get_model_evaluation_file_name()){
    path = get_model_eval_file_path()
    file_name = paste0(path, '/model_evaluation_scatter.pdf')
    eval_data = fread(file_path)
    eval_data[, length := motif_length_5_end + motif_length_3_end]

    plot = ggplot(eval) +
        geom_point(aes(x = length, y = log_loss, color = motif_length_5_end), size = 5) +
        theme_cowplot(font_family = 'Arial') + 
        xlab('Total motif length') +
        ylab('Conditional log loss') +
        theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
        guides(fill = guide_colourbar(barheight = 14)) + 
        scale_color_viridis_c(name = '5\' motif length')+
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 

    ggsave(file_name, plot = plot, width = 10, height = 7, units = 'in', dpi = 750, device = cairo_pdf)
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
