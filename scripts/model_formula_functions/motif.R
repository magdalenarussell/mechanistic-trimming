stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA

get_start_list <- function(motif_data, trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)
    all_pos = c()
    for (i in seq(length(trims))){
        all_pos = c(all_pos, get_positions(trims[i]))
    }
    positions = length(all_pos)
    start_list = rep(0, positions*3)
    return(start_list)
}

get_model_formula <- function(trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    motif_positions = c()
    for (i in seq(length(trims))){
        motif_positions = c(motif_positions, get_positions(trims[i]))
    }

    motif_positions_together = paste(motif_positions, collapse = ' + ')

    gene_groups = paste(paste0(genes, '_group'), collapse = ' ,')

    formula = formula(paste0('cbind(weighted_observation, interaction(', gene_groups, ', subject)) ~ ', motif_positions_together))
    return(formula)
}

process_single_data_for_model_fit <- function(group_motif_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    return(group_motif_data)
}
