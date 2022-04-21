stopifnot(LEFT_NUC_MOTIF_COUNT == 0 & RIGHT_NUC_MOTIF_COUNT == 0)

DATA_GROUP <<- 'ungrouped'

get_start_list <- function(motif_data){
    start_list = NULL
    return(start_list)
}

source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_base_count.R'))

get_model_formula <- function(){
    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ as.numeric(trim_length) + right_gc_prop'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    together = process_for_two_side_base_count(group_motif_data, left_nuc_count = 0)
    together[trim_length > 2, right_gc_prop := right_base_count_GC/(right_base_count_AT + right_base_count_GC)]
    together[trim_length <= 2, right_gc_prop := 0]
    return(together)
}
