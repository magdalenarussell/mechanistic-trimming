get_predicted_dist_figure_file_name <- function(gene_name, subgroup = NULL){
    new_gene_name = str_replace(gene_name, '/', '_')
    filename = paste0(new_gene_name, '_dist.pdf')
    return(filename)
}

get_coef_heatmap_file_name <- function(subgroup = NULL){
    filename = paste0('heatmap.pdf')
    return(filename)
}

get_residual_figure_file_name <- function(gene_name, subgroup = NULL){
    new_gene_name = str_replace(gene_name, '/', '_')
    filename = paste0(new_gene_name, '_residuals.pdf')
    return(filename)
}

get_all_residual_figure_file_name <- function(subgroup = NULL){
    filename = paste0('residuals.pdf')
    return(filename)
}

plot_predicted_trimming_dists <- function(data, gene_name){
    file_path = get_predicted_dist_figure_file_path()
    file_name = get_predicted_dist_figure_file_name(gene_name)
    complete_path = file.path(file_path, file_name)
    plot_predicted_trimming_dists_single_group(data, gene_name, complete_path)
}

plot_model_coefficient_heatmap <- function(model_coef_matrix, with_values = FALSE){
    file_path = get_coef_heatmap_file_path()
    file_name = get_coef_heatmap_file_name()
    complete_path = file.path(file_path, file_name)
    plot_model_coefficient_heatmap_single_group(model_coef_matrix, complete_path, with_values)
}

plot_model_residual_boxplot <- function(data, gene_name){
    file_path = get_residual_figure_file_path()
    file_name = get_residual_figure_file_name(gene_name)
    complete_path = file.path(file_path, file_name)
    plot_model_residual_boxplot_single_group(data, gene_name, complete_path)
}

plot_all_model_residuals <- function(data){
    file_path = get_all_residual_figure_file_path()
    file_name = get_all_residual_figure_file_name()
    complete_path = file.path(file_path, file_name)
    plot_all_model_residuals_plot(data, complete_path)
}
