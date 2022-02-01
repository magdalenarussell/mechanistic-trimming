stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)

DATA_GROUP <<- 'ungrouped'

get_start_list <- function(motif_data){
    dists = length(unique(motif_data$trim_length)) - 1
    positions = length(get_positions())
    start_list = rep(0, dists + positions*3 + 1)
    return(start_list)
}

source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/terminal_gc_content.R'))

get_model_formula <- function(){
    motif_positions = get_positions() 
    motif_positions_together = paste(motif_positions, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', motif_positions_together, ' + as.factor(trim_length) + terminal_gc_content'))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    together = process_for_terminal_gc_content(group_motif_data)
    return(together)
}
