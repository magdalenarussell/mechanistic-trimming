get_feature <- function(predicted_trims){
    p_gene_marg = unique(predicted_trims[, c('gene', 'p_gene')])
    setnames(p_gene_marg, 'p_gene', 'feature')
    return(p_gene_marg)
}
