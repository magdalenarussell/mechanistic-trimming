stopifnot((LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0) & (LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT) < 3)

DATA_GROUP <<- 'ungrouped'

get_start_list <- function(motif_data){
    unique_motifs = length(unique(motif_data$motif))
    dists = length(unique(motif_data$trim_length)) - 1
    start_list = rep(0, dists + unique_motifs-1)
    return(start_list)
}

get_model_formula <- function(){
    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ as.factor(motif) + as.factor(trim_length)'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    return(group_motif_data)
}
