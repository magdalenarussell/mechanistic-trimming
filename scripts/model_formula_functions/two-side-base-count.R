stopifnot(RIGHT_NUC_MOTIF_COUNT == 0 & LEFT_NUC_MOTIF_COUNT == 0)

DATA_GROUP <<- 'ungrouped'

get_start_list <- function(motif_data){
    start_list = NULL 
    return(start_list)
}

source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_base_count.R'))

get_model_formula <- function(){
    left_vars = get_all_base_variables('left')
    right_vars = get_all_base_variables('right')

    left_vars_collapse = paste(left_vars, collapse = ' + ')
    right_vars_collapse = paste(right_vars, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', left_vars_collapse, ' + ', right_vars_collapse))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    together = process_for_two_side_base_count(group_motif_data)
    return(together)
}
