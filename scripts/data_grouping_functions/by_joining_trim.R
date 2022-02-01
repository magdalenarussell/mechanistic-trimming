condense_tcr_data <- function(tcr_dataframe){
    cols = c('trimmed_seq', 'sequence', 'gene', paste(TRIM_TYPE), paste0(GENE_NAME), paste0(JOINING_TRIM))
    tcr_dataframe = tcr_dataframe[,..cols]
    setnames(tcr_dataframe, 'sequence', 'whole_seq')
    setnames(tcr_dataframe, TRIM_TYPE, 'trim_length')
    #condense data by gene, trim, etc.
    condensed_tcr = tcr_dataframe[, .N, by = .(trimmed_seq, whole_seq, gene, trim_length, get(GENE_NAME), get(JOINING_TRIM))]
    setnames(condensed_tcr, 'N', 'count')
    setnames(condensed_tcr, 'get', GENE_NAME)
    setnames(condensed_tcr, 'get.1', JOINING_TRIM)
    return(condensed_tcr)
}

sum_trim_observations <- function(condensed_tcr_dataframe){
    summed = condensed_tcr_dataframe[, sum(count), by = .(gene, get(JOINING_TRIM), trim_length, left_nucs, right_nucs)]
    setnames(summed, 'V1', 'count')
    setnames(summed, 'get', 'joining_trim')
    return(summed)
}

