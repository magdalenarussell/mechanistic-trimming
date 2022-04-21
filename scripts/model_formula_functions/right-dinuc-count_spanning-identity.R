stopifnot(RIGHT_NUC_MOTIF_COUNT == 0 & LEFT_NUC_MOTIF_COUNT == 0)

DATA_GROUP <<- 'ungrouped'

get_start_list <- function(motif_data){
    start_list = NULL 
    return(start_list)
}

source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_dinuc_count.R'))

get_model_formula <- function(){
    right_vars = get_all_dinucs_variables('right')

    right_vars_collapse = paste(right_vars, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', right_vars_collapse, ' + as.factor(span_identity)'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    together = process_for_two_side_dinuc_count(group_motif_data, left_nuc_count = 0)
    return(together)
}
