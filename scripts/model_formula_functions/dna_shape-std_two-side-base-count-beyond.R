stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)
stopifnot(RIGHT_NUC_MOTIF_COUNT < 3)

DATA_GROUP <<- 'ungrouped'

get_start_list <- function(motif_data, shapes = c('MGW', 'HelT', 'Roll', 'EP', 'ProT')){
    start_list = NULL
    return(start_list)
}

source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/dna_shape.R'))
source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_base_count.R'))

get_model_formula <- function(shapes = c('MGW', 'HelT', 'Roll', 'EP', 'ProT')){
    all_positions = c()
    for (shape in shapes){
        shape_positions = get_dna_shape_names_positions(shape)
        all_positions = c(all_positions, shape_positions)
    }

    all_positions = paste0(all_positions, '_std')
    shape_positions_together = paste(all_positions, collapse = ' + ')

    left_vars = get_all_base_variables('left')
    right_vars = get_all_base_variables('right')

    left_vars_collapse = paste(left_vars, collapse = ' + ')
    right_vars_collapse = paste(right_vars, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', right_vars_collapse, ' + ', left_vars_collapse, ' + ', shape_positions_together))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    together = process_for_two_side_base_count(group_motif_data, beyond_motif = TRUE)
    together = process_for_dna_structure(together, standardize = TRUE)
    stopifnot(nrow(together) == nrow(group_motif_data))
    return(together)
}
