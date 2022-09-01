stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)

DATA_GROUP <<- 'ungrouped'

get_start_list <- function(motif_data){
    positions = length(get_positions())
    start_list = rep(0, positions*3)
    return(start_list)
}

get_model_formula <- function(){
    motif_positions = get_positions() 
    motif_positions_together = paste(motif_positions, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', motif_positions_together))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data, whole_nucseq = get_oriented_whole_nucseqs()){
    return(group_motif_data)
}
