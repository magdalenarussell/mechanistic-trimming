calculate_subject_gene_weight <- function(compiled_data){
    compiled_data[, total_tcr := sum(count), by = subject]
    compiled_data[, gene_count := sum(count), by = .(gene, subject)]
    compiled_data$p_gene = compiled_data$gene_count/compiled_data$total_tcr
    compiled_data$p_trim_given_gene = compiled_data$count/compiled_data$gene_count
    p_subject = 1/length(unique(compiled_data$subject))
    compiled_data$weighted_observation = compiled_data$p_trim_given_gene*compiled_data$p_gene*p_subject
    compiled_data$gene_weight_type = 'p_gene_given_subject' 
    return(compiled_data)
}
