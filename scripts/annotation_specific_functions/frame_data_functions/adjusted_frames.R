source(paste0(MOD_PROJECT_PATH,'/scripts/annotation_specific_functions/frame_data_functions/all_frames.R'))

TRIMMING_LIGATION_DOMAIN <<- 'adjusted_mh'

get_frames_data <- function(){
    # get all possible configurations
    all_frames = get_all_frames_data()

    # get possible configs after MH adjustment
    zero_configs = all_frames[ligation_mh == 0]

    # adjust configs for MH
    converted_configs = get_possible_mh(zero_configs, keep_gene_seqs = FALSE)
    converted_configs = count_mh_bordering_trim(converted_configs)
    converted_configs = reassign_trimming_sites_with_mh(converted_configs)

    cols = c('v_gene', 'j_gene', 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'adjusted_v_trim', 'adjusted_j_trim', 'ligation_mh')

    ## get zero/nonzero MH total possible annotations after MH adjustment
    converted_condensed = converted_configs[, .N, by = cols]
    setnames(converted_condensed, 'adjusted_v_trim', 'v_trim')
    setnames(converted_condensed, 'adjusted_j_trim', 'j_trim')

    # get full sequences
    converted_condensed = get_oriented_full_sequences(converted_condensed)

    converted_condensed = get_stop_codon_positions(converted_condensed)

    # Frame calculations
    ## j_frame is subtracted because the overall sequence length should have a count of j_frame excess nucleotides 
    converted_condensed[, overall_frame := (v_seq_len + j_seq_len - (v_trim + j_trim) - ligation_mh - j_frame) %% 3]
    converted_condensed[, frame_type := fifelse(overall_frame == 0, 'In', 'Out')]

    # Look for stop codons based on frame
    converted_condensed = get_stop_positions_with_frame(converted_condensed)

    # Subset data by columns
    cols3 = c('v_gene', 'j_gene', 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'v_trim', 'j_trim', 'ligation_mh', 'overall_frame', 'frame_type', 'frame_stop')
    final = converted_condensed[, ..cols3]
    return(final)
}
