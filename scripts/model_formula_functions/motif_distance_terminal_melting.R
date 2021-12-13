source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/terminal_melting.R'))

get_model_formula <- function(){
    motif_positions = get_positions() 
    motif_positions_together = paste(motif_positions, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', motif_positions_together, ' + as.factor(trim_length) + terminal_melting'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    together = process_for_terminal_melting(group_motif_data, 'simple')
    return(together)
}
