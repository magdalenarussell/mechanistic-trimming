get_model_formula <- function(){
    motif_positions = get_positions() 
    motif_positions_together = paste(motif_positions, collapse = ' + ')

    interactions = ''
    for (index in seq(length(motif_positions)-1)){
        dinuc = paste(c(motif_positions[index], motif_positions[index + 1]), collapse =  ' * ')
        if (nchar(interactions) > 1){
            interactions = paste(c(interactions, dinuc), collapse = ' + ')
        } else {
            interactions = dinuc
        }
    }

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', interactions))
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    return(group_motif_data)
}
