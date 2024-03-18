TRIMMING_LIGATION_DOMAIN <<- 'all_mh'

get_all_frames_data <- function(){
    # This function will only return frame data for annotations that are observable, after re-annotation to accommodate cases of maximal ligation-mh
    # Read only necessary columns and apply filter
    frames = fread('https://raw.githubusercontent.com/phbradley/conga/master/conga/tcrdist/db/combo_xcr.tsv', select = c('id', 'region', 'nucseq', 'organism', 'chain', 'frame'))[organism == 'human' & chain == substring(CHAIN_TYPE, 3, 3)]

    # Combine operations to reduce redundancy
    v = frames[region == 'V', .(v_gene = id, v_frame = frame, v_seq = nucseq)]
    j = frames[region == 'J', .(j_gene = id, j_frame = frame, j_seq = nucseq)]

    # Merge operations
    v$dummy = 1
    j$dummy = 1
    gene_pairs = merge(v, j, by = 'dummy', allow.cartesian = TRUE)

    # Get all trimming sites
    trims = data.table(expand.grid(v_trim = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND), 
                                   j_trim = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)))
    trims$dummy = 1

    # Merge genes and trims
    all = merge(gene_pairs, trims, by = 'dummy', allow.cartesian = TRUE)[, -c('dummy')]
    # Get sequence lengths and subset data to necessary columns
    all[, c('v_seq_len', 'j_seq_len') := .(nchar(v_seq), nchar(j_seq))]

    cols = c('v_gene', 'j_gene', 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'v_trim', 'j_trim', 'v_gene_sequence', 'j_gene_sequence')
    all = get_oriented_full_sequences(all)
    all = unique(all[, ..cols])

    # get possible ligation mh
    lig_mat = matrix(0, nrow = nrow(all), ncol = length(seq(1, 15)))

    for (overlap in seq(ncol(lig_mat))){
        lig = get_possible_ligation_mh_fixed_trim(all, overlap_count = overlap)
        lig_mat[, overlap] = lig
    }

    # get unique vals 
    unique_mh_list = apply(lig_mat, 1, unique)

    # combine
    adjusted_all = data.table(all, ligation_mh = unique_mh_list)

    # Expand the rows
    adjusted_grouped = as.data.table(adjusted_all[, unnest(.SD, cols = c("ligation_mh"))])

    # Get oriented full sequences and group genes by common features, also subset columns again
    cols2 = c(cols, 'ligation_mh')
    adjusted_grouped = unique(adjusted_grouped[, ..cols2])
    
    # Filter by trimming length
    adjusted_grouped = adjusted_grouped[v_trim <= UPPER_TRIM_BOUND & j_trim <= UPPER_TRIM_BOUND]

    # Get stop positions
    adjusted_grouped = get_stop_codon_positions(adjusted_grouped)

    # Frame calculations
    ## j_frame is subtracted because the overall sequence length should have a count of j_frame excess nucleotides 
    adjusted_grouped[, overall_frame := (v_seq_len + j_seq_len - (v_trim + j_trim) - ligation_mh - j_frame) %% 3]
    adjusted_grouped[, frame_type := fifelse(overall_frame == 0, 'In', 'Out')]

    # Look for stop codons based on frame
    adjusted_grouped = get_stop_positions_with_frame(adjusted_grouped)

    # Subset data by columns
    cols3 = c('v_gene', 'j_gene', 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'v_trim', 'j_trim', 'ligation_mh', 'overall_frame', 'frame_type', 'frame_stop')
    final = adjusted_grouped[, ..cols3]
    return(final)
}

get_frames_data <- function(){
    return(get_all_frames_data())
}
