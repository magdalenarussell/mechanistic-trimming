condense_tcr_data <- function(tcr_dataframe, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    cols = c(paste0(trim_type, 'med_seq'), paste0(gene_type, '_sequence'), paste0(gene_type, '_group'), paste(trim_type), paste0(JOINING_GENE))
    tcr_dataframe = tcr_dataframe[,..cols]
    setnames(tcr_dataframe, paste0(gene_type, '_sequence'), paste0(gene_type, '_whole_seq'))
    #condense data by gene, trim, etc.
    cols = c(paste0(trim_type, 'med_seq'), paste0(gene_type, '_whole_seq'), trim_type, paste0(gene_type, '_group'), JOINING_GENE)
    condensed_tcr = tcr_dataframe[, .N, by = cols]
    setnames(condensed_tcr, 'N', 'count')
    return(condensed_tcr)
}

sum_trim_observations <- function(condensed_tcr_dataframe, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    cols = c(paste0(gene_type, '_group'), JOINING_GENE, trim_type, paste0(trim_type, '_left_nucs'), paste0(trim_type, '_right_nucs'))
    summed = condensed_tcr_dataframe[, sum(count), by = cols]
    setnames(summed, 'V1', 'count')
    return(summed)
}

filter_motif_data_for_possible_sites <- function(motif_data){
    if (JOINING_INSERT == 'vj_insert'){
        frame_data = get_frames_data()
    }
    return(motif_data)
}
