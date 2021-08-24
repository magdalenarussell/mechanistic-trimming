UPPER_TRIM_BOUND <<- 18
#MR TODO change to 1 once we are done with model eval
LOWER_TRIM_BOUND <<- 2 
PNUC_COUNT <<- 0

get_motifs <- function(tcr_dataframe, subject_id){
    motif_data = general_get_motifs(tcr_dataframe, subject_id)
    return(motif_data)
}
 
set_contrasts <- function(group_motif_data, ref_base = 'A'){
    positions = get_positions()
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
    return(group_motif_data)
}



