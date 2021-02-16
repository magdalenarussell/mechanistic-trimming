#TODO adjust for bounded case
#Not sure what TODO for short v-genes here...some v genes are only 3nt long at
#full length, so if they are trimmed at all, we can't have a "bounded" motif to
#the left...
subset_data_by_trim_length_bound <- function(trim_length, tcr_dataframe, cdr3_variable, lower_bound = NA, upper_bound = NA){
    tcr_dataframe_subset = tcr_dataframe[get(TRIM_TYPE) == trim_length]
    temp_subset_columns = c(GENE_NAME, 'gene_full_length_seq', cdr3_variable, TRIM_TYPE)
    temp_subset = tcr_dataframe_subset[,.N,by = temp_subset_columns]
    if (!is.na(lower_bound)){
        temp_subset = temp_subset[nchar(get(cdr3_variable)) > lower_bound]
    } 
    if (!is.na(upper_bound)){
        temp_subset = temp_subset[nchar(get(cdr3_variable)) < upper_bound]
    }
    return(temp_subset)
}

create_PFM_by_trim_length <- function(trim_length, tcr_dataframe){
    cdr3_variable = paste0('cdr3_nucseq_from_', substring(TRIM_TYPE, 1, 1))
    temp_subset_full_motif = subset_data_by_trim_length_bound(trim_length, tcr_dataframe, cdr3_variable, lower_bound = 3)
    motif_vector = c()
    for (row in seq(1:nrow(temp_subset_full_motif))){
        motif = get_motif_context(whole_gene_seq = temp_subset_full_motif[row,]$gene_full_length_seq, 
                                  trimmed_gene_seq = temp_subset_full_motif[row,][[cdr3_variable]])
        motif_vector = c(motif_vector, rep(toString(motif), temp_subset_full_motif[row,]$N))
    }
    motif_vector_DNAString = DNAStringSet(motif_vector)
    PFM = consensusMatrix(motif_vector_DNAString,baseOnly=TRUE)
#TODO remove short motif section...update with '-' addition in motif
    #function...remove 'other' row in resulting PFM
    temp_subset_short_motif = subset_data_by_trim_length_bound(trim_length, tcr_dataframe, cdr3_variable, upper_bound = 4)
    #TODO test if this works for short vtrim=1 edge case
    for (row in seq(1:nrow(temp_subset_short_motif))){
        trimmed_gene_seq = temp_subset_short_motif[row,][[cdr3_variable]]
        #TODO for now, skip instances where v gene is fully trimmed back...
        if (trimmed_gene_seq == ''){
            next
        }
        motif = get_motif_context(whole_gene_seq = temp_subset_short_motif[row,]$gene_full_length_seq, 
                                  trimmed_gene_seq) 
        short_motif_vector = rep(toString(motif), temp_subset_short_motif[row,]$N)
        temp_short_motif_vector = DNAStringSet(short_motif_vector)
        
        temp_PFM = consensusMatrix(temp_short_motif_vector,baseOnly=TRUE)
        empty_position_matrix = matrix(0, ncol = LEFT_NUC_MOTIF_COUNT-nchar(trimmed_gene_seq), nrow = 5)
        complete_PFM = cbind(empty_position_matrix, temp_PFM)
        PFM = PFM + complete_PFM
    }
    return(PFM)
}

