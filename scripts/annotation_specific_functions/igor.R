#TODO update this
extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][3]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_")[[1]][1]
    return(localID)
}

general_get_motifs <- function(tcr_dataframe, subject_id){
    #filter data by trim bounds
    tcr_dataframe = tcr_dataframe[get(TRIM_TYPE) >= LOWER_TRIM_BOUND & get(TRIM_TYPE) <= UPPER_TRIM_BOUND]
    tcr_dataframe[, trimmed_seq := substring(sequence, 1, nchar(sequence) - get(TRIM_TYPE))]
    cols = c('trimmed_seq', 'sequence', 'gene', paste(TRIM_TYPE), paste0(GENE_NAME))
    tcr_dataframe = tcr_dataframe[,..cols]
    setnames(tcr_dataframe, 'sequence', 'whole_seq')
    setnames(tcr_dataframe, TRIM_TYPE, 'trim_length')
    #condense data by gene, trim, etc.
    condensed_tcr = tcr_dataframe[, .N, by = .(trimmed_seq, whole_seq, gene, trim_length, get(GENE_NAME))]
    setnames(condensed_tcr, 'N', 'count')
    setnames(condensed_tcr, 'get', GENE_NAME)
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

compile_motifs_for_subject <- function(file_path){
    temp_data = fread(file_path)
    if (GENE_NAME == 'd_gene'){
        temp_data = temp_data[d_gene != '-']
    }
    subject_id = extract_subject_ID(file_path)
    
    output_location = get_subject_motif_output_location() 
    dir.create(output_location, recursive = TRUE, showWarnings = FALSE)
    whole_nucseq = fread(get(paste0('WHOLE_NUCSEQS_', ANNOTATION_TYPE)))[, -c('name')]
    temp_data = merge(temp_data, whole_nucseq, by.x = GENE_NAME, by.y = 'gene')
    gene_seqs = whole_nucseq[substring(gene, 4, 4) == toupper(substring(GENE_NAME, 1, 1))]
    setnames(gene_seqs, 'gene', GENE_NAME)
    setnames(gene_seqs, 'sequence', 'sequences')
    gene_groups = get_common_genes_from_seqs(gene_seqs)
    together = merge(temp_data, gene_groups, by = GENE_NAME)

    motif_data = get_motifs(together, subject_id)

    fwrite(motif_data, file = file.path(output_location, paste0(subject_id, '.tsv')), sep = '\t')
}


