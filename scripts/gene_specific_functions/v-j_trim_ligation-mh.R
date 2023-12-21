stopifnot(LOCUS == 'alpha')
stopifnot(INSERTIONS == 'zero')

get_oriented_whole_nucseqs <- function(whole_nucseq = get_whole_nucseqs()){
    # reorient sequence so that it is 5 -> 3 on the actual trimmed strand
    whole_nucseq[substring(gene, 4, 4) == 'J', sequence := unlist(lapply(sequence, function(x) as.character(reverseComplement(DNAString(x)))))]
    whole_nucseq[substring(gene, 4, 4) == 'J', j_gene_sequence := sequence]
    whole_nucseq[substring(gene, 4, 4) == 'V', v_gene_sequence := sequence]
    return(whole_nucseq[, -c('sequence')])
}

source(paste0(MOD_PROJECT_PATH,'/scripts/gene_count_specific_functions/double.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/mh_functions.R'))

filter_motif_data_for_possible_sites <- function(motif_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME){
    # Ensure that insertions are set to 'zero'
    stopifnot(INSERTIONS == 'zero')

    # Define frame-related columns
    frame_cols = c('frame_type', 'frame_stop', 'overall_frame')

    # Remove frame-related columns if present in motif_data
    if (any(frame_cols %in% colnames(motif_data))){
        cols = colnames(motif_data)[!(colnames(motif_data) %in% frame_cols)]
        motif_data = motif_data[, ..cols]
    }

    # Retrieve gene order based on gene_type
    genes = get_gene_order(gene_type)

    # Read frame data and filter for possible sites
    frame_data = read_frames_data()
    possible_sites = frame_data[frame_type == 'Out' | frame_stop == TRUE] 

    # Define columns for filtering possible sites
    cols = c(paste0(genes, '_group'), 'frame_type', 'overall_frame', 'frame_stop', 'v_trim', 'j_trim', 'ligation_mh')
    possible_sites_subset = unique(possible_sites[, ..cols])

    # Merge motif data with possible sites subset
    cols2 = c(paste0(genes, '_group'), 'v_trim', 'j_trim', 'ligation_mh')
    tog = merge(motif_data, possible_sites_subset, by = cols2)
    return(tog)
}

adjust_trimming_sites_for_ligation_mh <- function(tcr_dataframe){
    mh_dt = get_possible_mh(tcr_dataframe, keep_gene_seqs = FALSE)
    mh_dt = count_mh_bordering_trim(mh_dt)
    mh_dt = reassign_trimming_sites_with_mh(mh_dt)
    mh_dt[, v_trim := adjusted_v_trim]
    mh_dt[, j_trim := adjusted_j_trim]
    return(mh_dt)
}
