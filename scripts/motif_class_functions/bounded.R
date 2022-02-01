LOWER_TRIM_BOUND <<- 2 
PNUC_COUNT <<- 2
 
# TODO change this back once we are done comparing models
# LOWER_TRIM_BOUND <<- max(RIGHT_NUC_MOTIF_COUNT - 2, 1)

get_all_nuc_contexts <- function(tcr_dataframe, subject_id){
    stopifnot(RIGHT_NUC_MOTIF_COUNT - 2 <= LOWER_TRIM_BOUND)
    motif_data = general_get_all_nuc_contexts(tcr_dataframe, subject_id)
    return(motif_data)
}

set_contrasts <- function(group_motif_data, ref_base = 'A'){
    positions = get_positions()
    for (position in positions){
        group_motif_data[[position]] = as.factor(group_motif_data[[position]])
        group_motif_data[[position]] = relevel(group_motif_data[[position]], ref_base)
        contrasts(group_motif_data[[position]]) = contr.sum(4)
    }

    if (grepl('distance', MODEL_TYPE, fixed = TRUE)){
        trim_count = UPPER_TRIM_BOUND - LOWER_TRIM_BOUND + 1
        group_motif_data$trim_length = as.factor(group_motif_data$trim_length)
        contrasts(group_motif_data$trim_length) = contr.sum(trim_count) 
    }

    return(group_motif_data)
}


