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

filter_motif_data_for_possible_sites <- function(motif_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME, trim_type = TRIM_TYPE){
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
    possible_sites_subset = get_all_possible_sites(gene_type)

    # Merge motif data with possible sites subset
    cols2 = c(paste0(genes), 'v_trim', 'j_trim', 'ligation_mh')
    tog = merge(motif_data, possible_sites_subset, by = cols2)

    # fill in remaining unobserved, but possible sites 
    ## Note: we are not including annotations that can be adjusted to a ligation-mh scenario; while these sites are possible, we are ignoring them since we have moved their counts to the ligation-mh scenario
    
    filled_tog = fill_in_missing_possible_sites(possible_sites_subset, tog, trim_type, gene_type)
    return(filled_tog)
}

get_all_possible_sites <- function(gene_type = GENE_NAME){
    # Retrieve gene order based on gene_type
    genes = get_gene_order(gene_type)

    # Read frame data and filter for possible sites
    frame_data = read_frames_data()
    possible_sites = frame_data[frame_type == 'Out' | frame_stop == TRUE] 

    # Define columns for filtering possible sites
    cols = c(paste0(genes), 'frame_type', 'overall_frame', 'frame_stop', 'v_trim', 'j_trim', 'ligation_mh')
    possible_sites_subset = unique(possible_sites[, ..cols])
    return(possible_sites_subset)
}

get_missing_possible_sites <- function(possible_sites, filtered_motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    genes = get_gene_order(gene_type)
    trims = get_trim_order(trim_type)

    # Get observed sets of gene pairs
    possible_sites$gene_pair = paste0(possible_sites$v_gene, '_', possible_sites$j_gene)
    filtered_motif_data$gene_pair = paste0(filtered_motif_data$v_gene, '_', filtered_motif_data$j_gene) 
    possible_sites_subset = possible_sites[gene_pair %in% unique(filtered_motif_data$gene_pair)]

    # get unobserved scenarios
    unobserved = possible_sites_subset[!filtered_motif_data, on = colnames(possible_sites_subset)]
    return(unobserved)
}

fill_in_missing_possible_sites <- function(possible_sites, filtered_motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    # get unobserved scenarios
    unobserved = get_missing_possible_sites(possible_sites, filtered_motif_data, trim_type, gene_type)

    if (nrow(unobserved) == 0){
        return(filtered_motif_data)
    } else {
        # get NT context for these scenarios
        for (i in seq(length(trims))){
            cols = c(paste0(genes[i]), trims[i], paste0(trims[i], '_left_nucs'), paste0(trims[i], '_right_nucs'))
            common_cols = c(paste0(genes[i]), trims[i])
            subset = unique(filtered_motif_data[, ..cols])
            unobserved = merge(unobserved, subset, by = common_cols, all.x = TRUE)
        }

        unobserved$subject = unique(filtered_motif_data$subject)
        unobserved$gene_type = unique(filtered_motif_data$gene_type)
        unobserved[[paste0(trim_type, '_observed')]] = FALSE
        unobserved$count = 0

        together = rbind(filtered_motif_data, unobserved, fill = TRUE)
        return(together)
    }
}

adjust_trimming_sites_for_ligation_mh <- function(tcr_dataframe){
    mh_dt = get_possible_mh(tcr_dataframe, keep_gene_seqs = FALSE)
    mh_dt = count_mh_bordering_trim(mh_dt)
    mh_dt = reassign_trimming_sites_with_mh(mh_dt)
    mh_dt[, v_trim := adjusted_v_trim]
    mh_dt[, j_trim := adjusted_j_trim]
    return(mh_dt)
}
