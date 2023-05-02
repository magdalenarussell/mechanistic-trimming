condense_tcr_data <- function(tcr_dataframe){
    cols = c('trimmed_seq', 'sequence', 'gene', paste(TRIM_TYPE), paste0(GENE_NAME), paste0(JOINING_GENE))
    tcr_dataframe = tcr_dataframe[,..cols]
    setnames(tcr_dataframe, 'sequence', 'whole_seq')
    setnames(tcr_dataframe, TRIM_TYPE, 'trim_length')
    #condense data by gene, trim, etc.
    condensed_tcr = tcr_dataframe[, .N, by = .(trimmed_seq, whole_seq, gene, trim_length, get(GENE_NAME), get(JOINING_GENE))]
    setnames(condensed_tcr, 'N', 'count')
    setnames(condensed_tcr, 'get', GENE_NAME)
    setnames(condensed_tcr, 'get.1', JOINING_GENE)
    return(condensed_tcr)
}

sum_trim_observations <- function(condensed_tcr_dataframe){
    summed = condensed_tcr_dataframe[, sum(count), by = .(gene, get(JOINING_GENE), trim_length, left_nucs, right_nucs)]
    setnames(summed, 'V1', 'count')
    setnames(summed, 'get', 'joining_gene')
    return(summed)
}

filter_motif_data_for_possible_sites <- function(motif_data){
    if (JOINING_INSERT == 'vj_insert'){
        frame_data = get_frames_data()
        inframe = frame_data[frame_type == 'In']
    }
    return(motif_data)
}
