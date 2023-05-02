stopifnot(LEFT_NUC_MOTIF_COUNT == 0 & RIGHT_NUC_MOTIF_COUNT == 0)

get_start_list <- function(motif_data){
    start_list = rep(0, 1)
    return(start_list)
}

get_model_formula <- function(){
    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ as.numeric(trim_length)'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    return(group_motif_data)
}


