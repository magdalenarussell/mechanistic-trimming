calculate_subject_gene_weight <- function(subject_data){
    total_tcr = nrow(subject_data[observed == TRUE])
    gene_counts = subject_data[observed == TRUE, .N, by = .(gene)]
    gene_counts$gene_weight = gene_counts$N/total_tcr
    together = merge(subject_data, gene_counts[, c('gene', 'gene_weight')], by = 'gene')
    return(together)
}


