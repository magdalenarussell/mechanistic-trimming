calculate_subject_gene_weight <- function(compiled_data){
    compiled_data[, total_tcr := sum(count), by = subject]
    compiled_data[, gene_count := sum(count), by = .(gene, subject)]

    gene_counts = compiled_data[observed == TRUE, sum(count), by = .(gene, subject, total_tcr)]
    gene_counts$p_gene = gene_counts$V1/gene_counts$total_tcr

    total_subjects = length(unique(gene_counts$subject))
    subject_sum = gene_counts[, sum(p_gene), by = gene]
    subject_sum$p_gene = subject_sum$V1/total_subjects
    subject_sum$gene_weight_type = 'p_gene_marginal' 
    together = merge(compiled_data[, -c('p_gene', 'gene_weight_type')], subject_sum[, c('gene', 'p_gene', 'gene_weight_type')], by = 'gene')
    together$p_trim_given_gene = together$count/together$gene_count
    p_subject = 1/total_subjects
    together$weighted_observation = together$p_trim_given_gene * together$p_gene * p_subject
    return(together)
}
