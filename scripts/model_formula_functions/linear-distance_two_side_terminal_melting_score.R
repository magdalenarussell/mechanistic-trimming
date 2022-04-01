stopifnot(RIGHT_NUC_MOTIF_COUNT == 0 & LEFT_NUC_MOTIF_COUNT == 0)

DATA_GROUP <<- 'ungrouped'

get_start_list <- function(motif_data){
    start_list = rep(0, 3)
    return(start_list)
}

source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_terminal_melting.R'))

get_model_formula <- function(){
    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ left_terminal_melting_score + right_terminal_melting_score + as.numeric(trim_length)'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    together = process_for_two_side_terminal_melting(group_motif_data, 'simple')
    together = transform_terminal_melting_to_score(together)
    return(together)
}
