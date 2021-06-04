get_predicted_dist_figure_file_name <- function(gene){
    path = file.path(PROJECT_PATH, 'plots', 'predicted_trimming_distributions')
    dir.create(path)

    complete_path = paste0(path, '/', gene, '_dist.pdf')
    return(complete_path)
}

get_actual_trimming_dist_by_gene_subject <- function(actual_data, gene_name){
    subset = actual_data[gene == gene_name]
    return(subset)
}

plot_predicted_trimming_dists <- function(actual_data, predicted_data, gene_name){
    actual_dist = get_actual_trimming_dist_by_gene_subject(actual_data, gene_name)
    actual_dist = actual_dist[order(subject, trim_length)]

    predicted_data = predicted_data[gene == gene_name]
    N = paste(unique(predicted_data$total_gene_count))
    title = paste0(gene_name, ', N = ', N)
    plot = ggplot() +
        geom_line(data = actual_dist[count_subject_gene <= 100], aes(x = trim_length, y = p_trim_given_gene, group = subject), size = 2, alpha = 0.7, color = 'grey') +
        geom_line(data = actual_dist[count_subject_gene > 100], aes(x = trim_length, y = p_trim_given_gene, group = subject), size = 2, alpha = 0.7) +
        geom_line(data = predicted_data, aes(x = trim_length, y = predicted_prob), size = 3, alpha = 0.7, color = 'blue') +
        ggtitle(title) +
        xlab('Number of trimmed nucleotides') +
        ylab('Probability') +
        theme_cowplot(font_family = 'Arial') + 
        theme(legend.position = "none", text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 
    new_gene_name = str_replace(gene_name, '/', '_')
    plot_name = get_predicted_dist_figure_file_name(new_gene_name) 
    ggsave(plot_name, plot = plot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
    print(paste0('finished plot for ', gene_name))
    return(plot)
}
    

