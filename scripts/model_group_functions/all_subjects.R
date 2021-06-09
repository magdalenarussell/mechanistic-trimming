GENE_WEIGHT_TYPE <<- 'p_gene_given_subject'

calculate_subject_gene_weight <- function(subject_data){
    total_tcr = nrow(subject_data[observed == TRUE])
    gene_counts = subject_data[observed == TRUE, .N, by = .(gene)]
    gene_counts$gene_weight = gene_counts$N/total_tcr
    gene_counts$gene_weight_type = GENE_WEIGHT_TYPE
    together = merge(subject_data, gene_counts[, c('gene', 'gene_weight')], by = 'gene')
    return(together)
}


