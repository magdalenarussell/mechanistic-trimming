stopifnot(LEFT_NUC_MOTIF_COUNT == 0 & RIGHT_NUC_MOTIF_COUNT == 0)

stopifnot(DATA_GROUP == 'by_joining_gene')

get_start_list <- function(motif_data){
    genes = length(unique(motif_data$joining_gene)) - 1
    start_list = rep(0, genes)
    return(start_list)
}

get_model_formula <- function(){
    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ as.factor(joining_gene)')) 
    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    return(group_motif_data)
}
