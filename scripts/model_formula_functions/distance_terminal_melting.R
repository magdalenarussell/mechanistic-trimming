stopifnot(RIGHT_NUC_MOTIF_COUNT == 0)

DATA_GROUP <<- 'ungrouped'

get_start_list <- function(motif_data){
    dists = length(unique(motif_data$trim_length)) - 1
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
