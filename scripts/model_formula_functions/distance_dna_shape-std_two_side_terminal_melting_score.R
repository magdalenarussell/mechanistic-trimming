stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)

DATA_GROUP <<- 'ungrouped'

get_start_list <- function(motif_data, shapes = c('MGW', 'HelT', 'Roll', 'EP', 'ProT')){
    all_positions = c()
    for (shape in shapes){
        shape_positions = get_dna_shape_names_positions(shape)
        all_positions = c(all_positions, shape_positions)
    }
    positions = length(all_positions)
    dists = length(unique(motif_data$trim_length)) - 1

    start_list = rep(0, positions + 2 + dists)
    return(start_list)
}

source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/dna_shape.R'))
source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_terminal_melting.R'))

get_model_formula <- function(shapes = c('MGW', 'HelT', 'Roll', 'EP', 'ProT')){
    all_positions = c()
    for (shape in shapes){
        shape_positions = get_dna_shape_names_positions(shape)
        all_positions = c(all_positions, shape_positions)
    }

    all_positions = paste0(all_positions, '_std')
    shape_positions_together = paste(all_positions, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', shape_positions_together, ' + right_terminal_melting_score + left_terminal_melting_score + as.factor(trim_length)'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    together = process_for_dna_structure(group_motif_data, standardize = TRUE)
    together = process_for_two_side_terminal_melting(together, 'simple')
    together = transform_terminal_melting_to_score(together)
    stopifnot(nrow(together) == nrow(group_motif_data))
    return(together)
}
