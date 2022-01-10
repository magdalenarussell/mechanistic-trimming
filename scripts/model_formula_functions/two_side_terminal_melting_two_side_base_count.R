source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_terminal_melting.R'))
source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_base_count.R'))

get_model_formula <- function(){
    left_vars = get_all_base_variables('left')
    right_vars = get_all_base_variables('right')

    left_vars_collapse = paste(left_vars, collapse = ' + ')
    right_vars_collapse = paste(right_vars, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', left_vars_collapse, ' + ', right_vars_collapse, ' + left_terminal_melting + right_terminal_melting'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    together = process_for_two_side_base_count(group_motif_data)
    together = process_for_two_side_terminal_melting(together, 'simple')
    return(together)
}
