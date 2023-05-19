stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA

get_start_list <- function(motif_data, trim_type = TRIM_TYPE){
    positions = length(get_positions(trim_type))
    start_list = rep(0, positions*3)
    return(start_list)
}

get_model_formula <- function(trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    motif_positions = get_positions(trim_type) 
    motif_positions_together = paste(motif_positions, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(', gene_type, '_group, subject)) ~ ', motif_positions_together))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    return(group_motif_data)
}
