get_predicted_dist_figure_file_name <- function(gene_name, subgroup){
    new_gene_name = str_replace(gene_name, '/', '_')
    filename = paste0(subgroup, '_', new_gene_name, '_dist.pdf')
    return(filename)
}

get_coef_heatmap_file_name <- function(subgroup){
    filename = paste0(subgroup, '_heatmap.pdf')
    return(filename)
}

get_residual_figure_file_name <- function(gene_name, subgroup){
    new_gene_name = str_replace(gene_name, '/', '_')
    filename = paste0(subgroup, '_', new_gene_name, '_residuals.pdf')
    return(filename)
}

plot_predicted_trimming_dists <- function(data, gene_name){
    file_path = get_predicted_dist_figure_file_path()
    for (indiv in unique(data$subject)){
        subset_data = data[subject == indiv]
        file_name = get_predicted_dist_figure_file_name(gene_name, indiv)
        complete_path = file.path(file_path, file_name)
        plot_predicted_trimming_dists_single_group(subset_data, gene_name, complete_path)
        print(paste0('finished plotting for ', indiv))
    }
}

plot_model_coefficient_heatmap <- function(model_coef_matrix, with_values = FALSE){
    file_path = get_coef_heatmap_file_path()
    limits = range(sapply(model_coef_matrix[,-c('base', 'model_group')], range, na.rm = TRUE))/log(10)
    for (indiv in unique(model_coef_matrix$model_group)){
        subset_data = model_coef_matrix[model_group == indiv]
        file_name = get_coef_heatmap_file_name(indiv)
        complete_path = file.path(file_path, file_name)
        plot_model_coefficient_heatmap_single_group(subset_data, complete_path, with_values, limits)
        print(paste0('finished plotting for ', indiv))
    }
}

plot_model_residual_boxplot <- function(data, gene_name){
    file_path = get_residual_figure_file_path()
    for (indiv in unique(data$subject)){
        subset_data = data[subject == indiv]
        file_name = get_residual_figure_file_name(gene_name, indiv)
        complete_path = file.path(file_path, file_name)
        plot_model_residual_boxplot_single_group(subset_data, gene_name, complete_path)
        print(paste0('finished plotting for ', indiv))
    }
}
