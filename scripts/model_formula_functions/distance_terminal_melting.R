stopifnot(RIGHT_NUC_MOTIF_COUNT == 0)

get_start_list <- function(){
    dists = (UPPER_TRIM_BOUND - LOWER_TRIM_BOUND)
    start_list = rep(0, dists + 1)
    return(start_list)
}

source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/terminal_melting.R'))

get_model_formula <- function(){
    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ terminal_melting + as.factor(trim_length)'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    together = process_for_terminal_melting(group_motif_data, 'simple')
    return(together)
}
