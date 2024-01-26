get_pnucs <- function(whole_gene_nucseq, orient, pnuc_count){
    whole_gene_nucseq = DNAStringSet(whole_gene_nucseq)
    stopifnot(orient %in% c('top', 'bottom'))
    # top means top strand oriented 5 > 3
    # bottom means bottom strand oriented 3 > 5
    if (orient == 'top'){
        possible_pnucs = substring(as.character(reverseComplement(whole_gene_nucseq)),1, pnuc_count) 
    } else {
        possible_pnucs = substring(as.character(reverseComplement(whole_gene_nucseq)), nchar(whole_gene_nucseq) - pnuc_count + 1, nchar(whole_gene_nucseq))
    }
    return(possible_pnucs)
}

get_overlapping_regions_ligation_mh <- function(v_gene_top_seq, j_gene_bottom_seq, v_trim, j_trim, pnucs = 2){
    require(Biostrings)
    v_pnucs = get_pnucs(v_gene_top_seq, 'top', pnucs)    
    j_pnucs = get_pnucs(j_gene_bottom_seq, 'bottom', pnucs)    
    
    v_gene_top_seq_p = paste0(v_gene_top_seq, v_pnucs)
    j_gene_bottom_seq_p = paste0(j_pnucs, j_gene_bottom_seq)
    
    # get overlapping seqs when aligning trim sites
    v_overlap = substring(v_gene_top_seq_p, nchar(v_gene_top_seq_p) - (2*pnucs + v_trim + j_trim) + 1 , nchar(v_gene_top_seq_p))
    j_overlap = substring(j_gene_bottom_seq_p, 1, (2*pnucs + v_trim + j_trim))

    # get v_gene.j_trimmed overlaps
    vg_jt = substring(v_overlap, 1, pnucs+j_trim)
    jt_vg = substring(j_overlap, 1, pnucs+j_trim)

    # get j_gene.v_trimmed overlaps
    jg_vt = substring(j_overlap, pnucs+j_trim + 1)
    vt_jg = substring(v_overlap, pnucs+j_trim + 1)

    return(data.table('v_gene.j_trimmed' = vg_jt, 'j_trimmed.v_gene' = jt_vg, 'j_gene.v_trimmed' = jg_vt, 'v_trimmed.j_gene' = vt_jg))
}

get_mh_ligation_mh <- function(seq1, seq2, aligning_trim){
    stopifnot(aligning_trim %in% c('j_trim', 'v_trim'))
    stopifnot(all(nchar(seq1) == nchar(seq2)))
    
    compl = c('A' = 'T', 'T' = 'A', 'G' = 'C', 'C' = 'G')
    match = c('A' = 'A/T', 'T' = 'A/T', 'C' = 'G/C', 'G' = 'G/C')
    
    max_len = nchar(seq1)
    
    if (max_len == 0 | nchar(seq2) == 0){
        mh = data.table()
    } else {
        names = get_mh_colnames(aligning_trim, max_len)
        
        # Precompute constant values
        seq1_compl = compl[as.character(strsplit(seq1, NULL)[[1]])]
        seq1_match = match[as.character(strsplit(seq1, NULL)[[1]])]
        
        # Use vectorized operations to create the matrix
        mh_matrix = matrix('-', ncol = max_len)
        mh_matrix[, seq1_compl == strsplit(seq2, NULL)[[1]]] = seq1_match[seq1_compl == strsplit(seq2, NULL)[[1]]]
        
        # Convert the matrix to data.table
        mh = as.data.table(mh_matrix)
        colnames(mh) = names
    }
    return(mh)
}

get_mh_colnames <- function(aligning_trim, positions){
    if (aligning_trim == 'j_trim'){
        names = paste0(aligning_trim, '_mh_pos_', seq(positions, 1)) 
    } else {
        names = paste0(aligning_trim, '_mh_pos_', seq(1, positions)) 
    }
    return(names)
}

fill_in_missing_mh_positions <- function(mh_dt, max_len, aligning_trim){
    stopifnot(max_len >= ncol(mh_dt))

    all_names = get_mh_colnames(aligning_trim, max_len)
    not_present = all_names[!(all_names %in% colnames(mh_dt))]
    if (length(not_present) > 0){
        additional = as.data.table(matrix(ncol = length(not_present))) 
        colnames(additional) = not_present
        mh_dt = cbind(mh_dt, additional)
    }
    setcolorder(mh_dt, all_names)
    return(mh_dt)
}

get_mh_and_fill <- function(seq1, seq2, aligning_trim, max_len){
    return(fill_in_missing_mh_positions(get_mh_ligation_mh(seq1, seq2, aligning_trim), max_len, aligning_trim)) 
}

get_mh_dataframe <- function(data, aligning_trim, aligning_gene){
    col1 = paste0(aligning_trim, 'med.', aligning_gene)
    col2 = paste0(aligning_gene, '.', aligning_trim, 'med')
    stopifnot(col1 %in% colnames(data))
    stopifnot(col2 %in% colnames(data))

    cols = c(col1, col2)
    subset = unique(data[, ..cols]) 
    max = max(nchar(data[[col1]]), nchar(data[[col2]]))
    # max = UPPER_TRIM_BOUND + PNUC_COUNT
    names = get_mh_colnames(aligning_trim, max)

    if (nrow(subset) == 1){
        subset[, paste0(names) := data.table(matrix(mapply(get_mh_and_fill, get(col1), get(col2), aligning_trim, max), nrow = 1))]
    } else {
        subset[, paste0(names) := as.data.table(apply(t(mapply(get_mh_and_fill, get(col1), get(col2), aligning_trim, max)), 2, unlist))]
    }

    data = merge(data, subset, by = cols)
    return(data)
}

get_possible_mh <- function(data, keep_gene_seqs = FALSE){
    # get full gene sequences
    genes = get_whole_nucseqs()
    data = merge(data, genes, by.x = 'v_gene', by.y = 'gene')
    setnames(data, 'sequence', 'v_gene_sequence')
    data = merge(data, genes, by.x = 'j_gene', by.y = 'gene')
    setnames(data, 'sequence', 'j_gene_sequence')

    # convert j sequence to be the bottom strand (oriented 3'->5')
    require(Biostrings)
    data[, j_gene_sequence := as.character(complement(DNAStringSet(j_gene_sequence)))]

    # get overlapping regions
    data[, c("v_gene.j_trimmed", "j_trimmed.v_gene", "j_gene.v_trimmed", "v_trimmed.j_gene") := get_overlapping_regions_ligation_mh(v_gene_sequence, j_gene_sequence, v_trim, j_trim)]

    # get MH
    data = get_mh_dataframe(data, aligning_trim = 'j_trim', aligning_gene = 'v_gene')
    data = get_mh_dataframe(data, aligning_trim = 'v_trim', aligning_gene = 'j_gene')

    if (keep_gene_seqs == FALSE) {
        remove = c('v_gene_sequence', 'j_gene_sequence')
        cols = colnames(data)[!(colnames(data) %in% remove)]
        data = data[, ..cols]
    }
    return(data)
}

count_mh_bordering_trim <- function(mh_data){
    mh_data[, bordering_mh_j_trim := 0]
    mh_data[, bordering_mh_v_trim := 0]
    mh_data[, bordering_mh_j_trim_nt := '']
    mh_data[, bordering_mh_v_trim_nt := '']
    
    # max = UPPER_TRIM_BOUND + PNUC_COUNT
    max = 10

    for (type in c('v_trim', 'j_trim')){
        baseline = 0
        for (pos in seq(1, max)){
            col = paste0(type, '_mh_pos_', pos)
            border_col = paste0('bordering_mh_', type)
            mh_data[get(col) != '-' & !is.na(get(col)) & get(border_col) == baseline, paste0(border_col) := get(border_col) + 1]
            if (type == 'j_trim'){
                mh_data[get(col) != '-' & !is.na(get(col)) & get(border_col) == baseline + 1, paste0(border_col, '_nt') := paste0(get(col), get(paste0(border_col, '_nt')))]
            } else {
                mh_data[get(col) != '-' & !is.na(get(col)) & get(border_col) == baseline + 1, paste0(border_col, '_nt') := paste0(get(paste0(border_col, '_nt')), get(col))]
            }
            baseline = baseline + 1
        }
    } 
    
    mh_data[, total_bordering_mh := bordering_mh_v_trim + bordering_mh_j_trim]
    mh_data[, total_bordering_mh_nt := paste0(bordering_mh_j_trim_nt, bordering_mh_v_trim_nt)]

    # this is a 1-based index from the end of the V-gene (not including pnucs)
    mh_data[total_bordering_mh > 0, bordering_mh_v_index_start := v_trim + bordering_mh_j_trim]
    mh_data[total_bordering_mh > 0, bordering_mh_v_index_end := v_trim - bordering_mh_v_trim + 1]
    return(mh_data)
}

reassign_trimming_sites_with_mh <- function(mh_data){
    mh_data[, adjusted_v_trim := v_trim - bordering_mh_v_trim]
    mh_data[, adjusted_j_trim := j_trim - bordering_mh_j_trim]
    mh_data[, ligation_mh := total_bordering_mh]
    return(mh_data)
}

get_mh <- function(seq1, seq2){
    seq2_comp = as.character(complement(DNAStringSet(seq2)))
    mh = mcmapply(function(X,Y) sum(str_count(X,Y)), strsplit(seq1, ''), strsplit(seq2_comp, ''))
    return(mh)
}

get_fully_contiguous_mh <- function(seq1, seq2){
    noncontig_mh = get_mh(seq1, seq2)
    contig_mh = ifelse(noncontig_mh == nchar(seq1), noncontig_mh, 0)
    return(contig_mh)
}

get_overlapping_regions <- function(v_gene_top_seq, j_gene_bottom_seq, v_trim, j_trim, overlap_count, pnucs = 2, positions = c('up', 'mid', 'down')){
    # j_gene_bottom_seq is the bottom strand oriented 5 > 3
    # this reorientation will orient the bottom strand 3 > 5 so that it is complementary to the V-gene top strand
    j_gene_bottom_seq = reorient_j_bottom_strand(j_gene_bottom_seq)

    require(Biostrings)
    v_pnucs = get_pnucs(v_gene_top_seq, 'top', pnucs)    
    j_pnucs = get_pnucs(j_gene_bottom_seq, 'bottom', pnucs)    
    
    v_gene_top_seq_p = paste0(v_gene_top_seq, v_pnucs)
    j_gene_bottom_seq_p = paste0(j_pnucs, j_gene_bottom_seq)
    
    # get overlapping seqs when aligning trim sites
    v_overlap = substring(v_gene_top_seq_p, nchar(v_gene_top_seq_p) - (2*pnucs + v_trim + j_trim + overlap_count) + 1 , nchar(v_gene_top_seq_p))
    j_overlap = substring(j_gene_bottom_seq_p, 1, (2*pnucs + v_trim + j_trim + overlap_count))

    # get v_gene.j_trimmed overlaps
    vg_jt = substring(v_overlap, 1, pnucs+j_trim)
    jt_vg = substring(j_overlap, 1, pnucs+j_trim)

    # get j_gene.v_trimmed overlaps
    jg_vt = substring(j_overlap, pnucs+j_trim + 1 + overlap_count)
    vt_jg = substring(v_overlap, pnucs+j_trim + 1 + overlap_count)

    # get overlaps
    if (overlap_count > 0){
        v_mid = substring(v_overlap, pnucs + j_trim + 1, pnucs + j_trim + overlap_count)
        j_mid = substring(j_overlap, pnucs + j_trim + 1, pnucs + j_trim + overlap_count)
    } else {
        v_mid = ""
        j_mid = ""
    }

    up = data.table(vg_jt, jt_vg)
    down = data.table(jg_vt, vt_jg)
    mid = data.table(v_mid, j_mid)

    overlaps = c()
    for (pos in positions){
        overlaps = cbind(overlaps, get(pos))
    }
    names = get_overlap_names(overlap_count, positions)
    setnames(overlaps, colnames(overlaps), names)

    return(overlaps)
}

get_possible_ligation_mh_fixed_trim <- function(data, overlap_count){
    # get overlapping regions
    names = get_overlap_names(overlap_count, 'mid')
    subset = data[, get_overlapping_regions(v_gene_sequence, j_gene_sequence, v_trim, j_trim, overlap_count, positions = 'mid')]
    colnames(subset) = names

    # get MH
    pos = 'mid'
    v_seq_col = paste0('v_gene_', pos, '_overlap_', overlap_count)
    j_seq_col = paste0('j_gene_', pos, '_overlap_', overlap_count)

    n = get_fully_contiguous_mh(subset[[v_seq_col]], subset[[j_seq_col]])
    return(n)
}



