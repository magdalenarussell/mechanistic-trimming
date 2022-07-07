get_predicted_dist_figure_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'predicted_trimming_distributions')
    dir.create(path, recursive = TRUE)
    return(path)
}

get_coef_heatmap_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }

    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'model_coefficient_heatmaps')
    dir.create(path, recursive = TRUE)
    return(path)
}

get_residual_figure_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }

    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'predicted_trimming_residuals')

    dir.create(path, recursive = TRUE)
    return(path)
}

get_coef_variations_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }

    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'model_coefficient_variations')

    dir.create(path, recursive = TRUE)
    return(path)
}

get_all_residual_figure_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }

    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'ALL_predicted_trimming_residuals')
    dir.create(path, recursive = TRUE)
    return(path)
}

get_base_composition_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }

    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'model_base_composition')
    dir.create(path, recursive = TRUE)
    return(path)
}

get_model_eval_file_path <- function(type){
    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, paste0('model_evaluation_', type), paste0(TRIM_TYPE, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
    dir.create(path, recursive = TRUE)
    return(path)
}

get_individual_comparison_file_path <- function(negative_control = FALSE, negative_control_type = NULL){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }

    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'individual_comparison')
    if (isTRUE(negative_control)){
        stopifnot(!is.null(negative_control_type))
        path = file.path(path, paste0('negative_control', negative_control_type))
    }
    dir.create(path, recursive = TRUE)
    return(path)
}

get_analysis_path <- function(){
    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, 'analysis')
    dir.create(path, recursive = TRUE)
    return(path)
}

get_manuscript_path <- function(){
    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, 'manuscript')
    dir.create(path, recursive = TRUE)
    return(path)
}
