condense_tcr_data <- function(tcr_dataframe){
    cols = c('trimmed_seq', 'sequence', 'gene', paste(TRIM_TYPE), paste0(GENE_NAME), paste(JOINING_INSERT))
    tcr_dataframe = tcr_dataframe[,..cols]
    tcr_dataframe = tcr_dataframe[get(JOINING_INSERT) == 0]
    setnames(tcr_dataframe, 'sequence', 'whole_seq')
    setnames(tcr_dataframe, TRIM_TYPE, 'trim_length')
    #condense data by gene, trim, etc.
    cols = c('trimmed_seq', 'whole_seq', 'gene', 'trim_length', GENE_NAME, JOINING_INSERT)
    condensed_tcr = tcr_dataframe[, .N, by = cols]
    setnames(condensed_tcr, 'N', 'count')
    return(condensed_tcr)
}

sum_trim_observations <- function(condensed_tcr_dataframe){
    summed = condensed_tcr_dataframe[, sum(count), by = .(gene, trim_length, left_nucs, right_nucs)]
    setnames(summed, 'V1', 'count')
    return(summed)
}

filter_motif_data_for_possible_sites <- function(motif_data){
    return(motif_data)
}
