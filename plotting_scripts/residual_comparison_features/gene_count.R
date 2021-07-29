get_feature <- function(predicted_trims){
    gene_count_sum = predicted_trims[, sum(count), by = gene]
    setnames(gene_count_sum, 'V1', 'feature')
    return(gene_count_sum)
}
