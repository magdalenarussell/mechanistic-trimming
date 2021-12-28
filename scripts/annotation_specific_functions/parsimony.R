extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][3]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_")[[1]][3]
    return(localID)
}

general_get_motifs <- function(tcr_dataframe, subject_id){
    #filter data by trim bounds
    tcr_dataframe = tcr_dataframe[get(TRIM_TYPE) >= LOWER_TRIM_BOUND & get(TRIM_TYPE) <= UPPER_TRIM_BOUND]
    cdr3_variable = paste0('cdr3_nucseq_from_', substring(TRIM_TYPE, 1, 1))
    cols = c(paste(cdr3_variable), 'sequence', 'gene', paste(TRIM_TYPE), paste0(GENE_NAME))
    tcr_dataframe = tcr_dataframe[,..cols]
    colnames(tcr_dataframe) = c('trimmed_seq', 'whole_seq', 'gene', 'trim_length', GENE_NAME)
    #condense data by gene, trim, etc.
    condensed_tcr = tcr_dataframe[, .N, by = .(trimmed_seq, whole_seq, gene, trim_length, get(GENE_NAME))]
    colnames(condensed_tcr) = c('trimmed_seq', 'whole_seq', 'gene', 'trim_length', GENE_NAME, 'count')
    #get motifs
    motif_dataframe = condensed_tcr[,motif:=get_motif_context(whole_seq, trimmed_seq, trim_length)]
    motif_dataframe$observed = TRUE
    unobserved = get_unobserved_motifs(motif_dataframe)
    together = rbind(motif_dataframe[,-c(1,2)], unobserved)
    recondensed = together[, sum(count), by = .(gene, trim_length, motif)]
    setnames(recondensed, 'V1', 'count')
    recondensed[count == 0, observed := FALSE]
    recondensed[count != 0, observed := TRUE]

    recondensed$gene_type = GENE_NAME
    recondensed$subject = subject_id
    return(recondensed)
}

get_whole_nucseqs <- function(){
    whole_nucseq = fread(get(paste0('WHOLE_NUCSEQS_', ANNOTATION_TYPE)))
    setnames(whole_nucseqs, 'gene_names', 'gene')
    setnames(whole_nucseqs, 'sequences', 'sequence')
    return(whole_nucseq)
}

