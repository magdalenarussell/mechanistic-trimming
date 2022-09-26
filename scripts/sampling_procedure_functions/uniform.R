calculate_subject_gene_weight <- function(compiled_data){
    compiled_data[, total_tcr := sum(count), by = subject]
    compiled_data[, gene_count := sum(count), by = .(gene, subject)]
    
    total_subjects = length(unique(compiled_data$subject))
    
    unique_genes = compiled_data[trim_length == 2, .N, by = .(subject)]
    compiled_data = merge(compiled_data, unique_genes, by = 'subject')
    compiled_data[, p_gene := 1/N]
    compiled_data$gene_weight_type = 'uniform' 
    compiled_data$p_trim_given_gene = compiled_data$count/compiled_data$gene_count
    p_subject = 1/total_subjects
    compiled_data$weighted_observation = compiled_data$p_trim_given_gene * compiled_data$p_gene * p_subject
    return(compiled_data[, -c('N')])
}
