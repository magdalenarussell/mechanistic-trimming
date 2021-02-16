subset_data_by_trim_length_bound <- function(trim_length, tcr_dataframe, cdr3_variable){
    tcr_dataframe_subset = tcr_dataframe[get(TRIM_TYPE) == trim_length]
    temp_subset_columns = c(GENE_NAME, 'sequences', cdr3_variable, TRIM_TYPE)
    temp_subset = tcr_dataframe_subset[,.N,by = temp_subset_columns]
    return(temp_subset)
}

create_PFM_by_trim_length <- function(trim_length, tcr_dataframe){
    cdr3_variable = paste0('cdr3_nucseq_from_', substring(TRIM_TYPE, 1, 1))
    temp_subset = subset_data_by_trim_length_bound(trim_length, tcr_dataframe, cdr3_variable)
    motif_vector = c()
    for (row in seq(1:nrow(temp_subset))){
        motif = get_motif_context(whole_gene_seq = temp_subset[row,]$sequences, 
                                  trimmed_gene_seq = temp_subset[row,][[cdr3_variable]], 
                                  trim_length)
        motif_vector = c(motif_vector, rep(toString(motif), temp_subset[row,]$N))
    }
    motif_vector_DNAString = DNAStringSet(motif_vector)
    PFM = consensusMatrix(motif_vector_DNAString,baseOnly=TRUE)
    return(PFM)
}


