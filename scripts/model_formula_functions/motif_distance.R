stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)

get_start_list <- function(){
    dists = (UPPER_TRIM_BOUND - LOWER_TRIM_BOUND)
    positions = length(get_positions())
    start_list = rep(0, dists + positions*3)
    return(start_list)
}

get_model_formula <- function(){
    motif_positions = get_positions() 
    motif_positions_together = paste(motif_positions, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', motif_positions_together, ' + as.factor(trim_length)'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    return(group_motif_data)
}
