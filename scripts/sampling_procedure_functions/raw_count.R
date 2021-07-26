calculate_subject_gene_weight <- function(compiled_data){
    compiled_data$weighted_observation = compiled_data$count 
    compiled_data$gene_weight_type = GENE_WEIGHT_TYPE
    return(compiled_data)
}
