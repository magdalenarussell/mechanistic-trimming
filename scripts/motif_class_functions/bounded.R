get_unobserved_motifs <- function(tcr_dataframe, bound){
    tcr_dataframe_observed = tcr_dataframe[,.N, by = .(gene, trim_length)]
    unobserved_df = data.frame()
    for (gene_name in unique(tcr_dataframe_observed$gene)){
        observed = unique(tcr_dataframe_observed[gene == paste(gene_name)]$trim_length)
        unobserved = setdiff(seq(bound,UPPER_TRIM_BOUND), observed)
        unobserved_df = rbind(unobserved_df, data.frame(gene = rep(gene_name, length(unobserved)), trim_length = unobserved))
    }
    cols = c('whole_seq', 'gene', GENE_NAME)
    tcr_dataframe = unique(tcr_dataframe[,..cols])
    together = as.data.table(merge(unobserved_df, tcr_dataframe, by = 'gene'))
    together = together[, motif:= get_motif_context_unobserved(whole_seq, trim_length)] 
    together$observed = FALSE
    together$count = 0
    return(together[,-c('whole_seq')])
}

get_motifs <- function(tcr_dataframe, subject_id){
    bound = LOWER_TRIM_BOUND
    #filter data by trim bounds
    tcr_dataframe = tcr_dataframe[get(TRIM_TYPE) >= bound & get(TRIM_TYPE) <= UPPER_TRIM_BOUND]
    cdr3_variable = paste0('cdr3_nucseq_from_', substring(TRIM_TYPE, 1, 1))
    cols = c(paste(cdr3_variable), 'sequences', 'gene', paste(TRIM_TYPE), paste0(GENE_NAME))
    tcr_dataframe = tcr_dataframe[,..cols]
    colnames(tcr_dataframe) = c('trimmed_seq', 'whole_seq', 'gene', 'trim_length', GENE_NAME)
    #condense data by gene, trim, etc.
    condensed_tcr = tcr_dataframe[, .N, by = .(trimmed_seq, whole_seq, gene, trim_length, get(GENE_NAME))]
    colnames(condensed_tcr) = c('trimmed_seq', 'whole_seq', 'gene', 'trim_length', GENE_NAME, 'count')
    #get motifs
    motif_dataframe = condensed_tcr[,motif:=get_motif_context(whole_seq, trimmed_seq, trim_length)]
    motif_dataframe$observed = TRUE
    unobserved = get_unobserved_motifs(motif_dataframe, bound)
    together = rbind(motif_dataframe[,-c(1,2)], unobserved)
    together$gene_type = GENE_NAME
    together$subject = subject_id
    return(together)
}

