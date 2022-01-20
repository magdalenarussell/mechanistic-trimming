stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)

source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_terminal_melting.R'))

get_model_formula <- function(){
    motif_positions = get_positions() 
    motif_positions_together = paste(motif_positions, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', motif_positions_together, ' + left_terminal_melting + right_terminal_melting'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    together = process_for_two_side_terminal_melting(group_motif_data, "simple")
    return(together)
}
