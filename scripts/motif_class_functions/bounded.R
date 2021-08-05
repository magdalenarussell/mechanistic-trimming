UPPER_TRIM_BOUND <<- 18
LOWER_TRIM_BOUND <<- 2
# TODO change this back once we are done comparing models
# LOWER_TRIM_BOUND <<- max(RIGHT_NUC_MOTIF_COUNT - 2, 1)

get_motifs <- function(tcr_dataframe, subject_id){
    stopifnot(RIGHT_NUC_MOTIF_COUNT - 2 <= LOWER_TRIM_BOUND)
    motif_data = general_get_motifs(tcr_dataframe, subject_id)
    return(motif_data)
}

set_contrasts <- function(group_motif_data, ref_base = 'A'){
    positions = get_positions()
    for (position in positions){
        group_motif_data[[position]] = as.factor(group_motif_data[[position]])
        group_motif_data[[position]] = relevel(group_motif_data[[position]], ref_base)
        contrasts(group_motif_data[[position]]) = contr.sum(4)
    }
    return(group_motif_data)
}


