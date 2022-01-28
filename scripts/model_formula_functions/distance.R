stopifnot(LEFT_NUC_MOTIF_COUNT == 0 & RIGHT_NUC_MOTIF_COUNT == 0)

get_start_list <- function(){
    dists = (UPPER_TRIM_BOUND - LOWER_TRIM_BOUND)
    start_list = rep(0, dists)
    return(start_list)
}

get_model_formula <- function(){
    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ as.factor(trim_length)'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    return(group_motif_data)
}


