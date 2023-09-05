PNUC_COUNT <<- -1 

get_all_nuc_contexts <- function(tcr_dataframe, subject_id, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    motif_data = general_get_all_nuc_contexts(tcr_dataframe, subject_id, gene_type = gene_type, trim_type = trim_type)
    return(motif_data)
}

single_set_contrasts <- function(group_motif_data, ref_base = 'A', trim_type = TRIM_TYPE){
    positions = get_positions(trim_type)
    for (position in positions){
        group_motif_data[[position]] = as.factor(group_motif_data[[position]])
        group_motif_data[[position]] = relevel(group_motif_data[[position]], ref_base)
        unique_bases = length(unique(group_motif_data[[position]]))
        contrasts(group_motif_data[[position]]) = contr.sum(unique_bases)
        if (unique_bases > 4){
            contrasts = contrasts(group_motif_data[[position]])
            non_base_matrix_pos = which(!(rownames(contrasts) %in% c('A', 'T', 'C', 'G')))
            stopifnot(non_base_matrix_pos < length(rownames(contrasts)))
            contrasts[non_base_matrix_pos,] = rep(0, length(contrasts[non_base_matrix_pos,]))
            contrasts[,non_base_matrix_pos] = rep(0, length(contrasts[,non_base_matrix_pos]))
            contrasts(group_motif_data[[position]]) = contrasts
        }
    }

    if (grepl('distance', MODEL_TYPE, fixed = TRUE) & !grepl('linear-distance', MODEL_TYPE, fixed = TRUE)){
        trim_count = UPPER_TRIM_BOUND - LOWER_TRIM_BOUND + 1
        group_motif_data[[trim_type]] = as.factor(group_motif_data[[trim_type]])
        contrasts(group_motif_data[[trim_type]]) = contr.sum(trim_count) 
    }

    return(group_motif_data)
}



