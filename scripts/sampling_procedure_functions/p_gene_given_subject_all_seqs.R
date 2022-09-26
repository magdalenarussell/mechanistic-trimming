calculate_subject_gene_weight <- function(compiled_data){
    if (PRODUCTIVITY != 'both'){
        productive = PRODUCTIVITY
        PRODUCTIVITY <<- 'both'
        other_prod_data = aggregate_validation_data(VALIDATION_DATA_DIR)
        PRODUCTIVITY <<- productive
    } else {
        other_prod_data = compiled_data
    }
    cols = colnames(other_prod_data)[!(colnames(other_prod_data) == "count")]
    other_prod_data = other_prod_data[, sum(count), by = cols]
    setnames(other_prod_data, 'V1', 'count')
    other_prod_data[, total_tcr := sum(count), by = subject]
    other_prod_data[, gene_count := sum(count), by = .(gene, subject)]

    gene_counts = other_prod_data[observed == TRUE, sum(count), by = .(gene, subject, total_tcr)]
    gene_counts$p_gene = gene_counts$V1/gene_counts$total_tcr
    total_subjects = length(unique(gene_counts$subject))
    gene_counts$gene_weight_type = 'p_gene_given_subject_all_seqs' 
    together = merge(compiled_data[, -c('p_gene', 'gene_weight_type')], gene_counts[, c('subject', 'gene', 'p_gene', 'gene_weight_type')], by = c('subject', 'gene'))
    weights = together[trim_length == 2, sum(p_gene), by = subject]
    together = merge(together, weights)
    together[, p_gene := p_gene/V1]
    together$p_trim_given_gene = together$count/together$gene_count
    p_subject = 1/total_subjects
    together$weighted_observation = together$p_trim_given_gene * together$p_gene * p_subject
    return(together[, -c('V1')])
}
