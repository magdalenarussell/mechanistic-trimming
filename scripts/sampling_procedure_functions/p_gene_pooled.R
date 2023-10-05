calculate_subject_gene_weight <- function(compiled_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    genes = get_gene_order(gene_type)
    trims = get_trim_order(trim_type)

    if ('subject' %in% colnames(compiled_data)){
        params = get_parameter_vector(trims, genes)
        cols = c(paste0(genes, '_group'), trims, params)
        compiled_data = compiled_data[, sum(count), by = cols]
        setnames(compiled_data, 'V1', 'count')
    }
    compiled_data[, total_tcr := sum(count)]
    col = c(paste0(genes, '_group'))
    compiled_data[, paste0(gene_type, '_count') := sum(count), by = col]
    compiled_data[[paste0('p_', gene_type)]] = compiled_data[[paste0(gene_type, '_count')]]/compiled_data$total_tcr
    compiled_data[[paste0('p_', trim_type, '_given_', gene_type)]] = compiled_data$count/compiled_data[[paste0(gene_type, '_count')]]
    compiled_data$weighted_observation = compiled_data[[paste0('p_', trim_type, '_given_', gene_type)]]*compiled_data[[paste0('p_', gene_type)]]
    compiled_data$gene_weight_type = 'p_gene_pooled' 
    return(compiled_data)
}
