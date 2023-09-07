get_pnucs <- function(whole_gene_nucseq, orient, pnuc_count){
    whole_gene_nucseq = DNAStringSet(whole_gene_nucseq)
    stopifnot(orient %in% c('top', 'bottom'))
    if (orient == 'top'){
        possible_pnucs = substring(as.character(reverseComplement(whole_gene_nucseq)),1, pnuc_count) 
    } else {
        possible_pnucs = substring(as.character(reverseComplement(whole_gene_nucseq)), nchar(whole_gene_nucseq) - pnuc_count + 1, nchar(whole_gene_nucseq))
    }
    return(possible_pnucs)
}

reorient_j_bottom_strand <- function(j_gene_5_3_bottom_strand){
    # return j_gene bottom strand (which was previously oriented 5' > 3' for consistency with the v_gene) to be 3' > 5'
    require(stringi)
    return(stri_reverse(j_gene_5_3_bottom_strand))
}

get_overlap_names <- function(overlap_count){
    names = c(paste0('v_gene_up_overlap_', overlap_count), 
              paste0('j_gene_up_overlap_', overlap_count), 
              paste0('j_gene_down_overlap_', overlap_count),
              paste0('v_gene_down_overlap_', overlap_count),
              paste0('v_gene_mid_overlap_', overlap_count), 
              paste0('j_gene_mid_overlap_', overlap_count))
    return(names)
}

get_all_mh_prop_variables <- function(overlap_vector){
    pos = c('up', 'mid', 'down')
    vars = c()
    for (o in overlap_vector){
        for (p in pos){
            if (o == 0 & p == 'mid'){
                next
            }
            var = paste0('mh_prop_', p, '_overlap_', o) 
            vars = c(vars, var)
        }
    }
    return(vars)
}

get_overlapping_regions <- function(v_gene_top_seq, j_gene_bottom_seq, v_trim, j_trim, overlap_count, pnucs = 2){
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

    overlaps = data.table(vg_jt, jt_vg, jg_vt, vt_jg, v_mid, j_mid)
    names = get_overlap_names(overlap_count)
    setnames(overlaps, colnames(overlaps), names)

    return(overlaps)
}

get_mh <- function(seq1, seq2){
    seq2_comp = as.character(complement(DNAStringSet(seq2)))
    mh = mcmapply(function(X,Y) sum(str_count(X,Y)), strsplit(seq1, ''), strsplit(seq2_comp, ''))
    return(mh)
}

get_mh_prop <- function(seq1, seq2){
    mh = get_mh(seq1, seq2)
    lengths = nchar(seq1) 
    mh_prop = mh/lengths
    return(mh_prop)
}

get_mh_prop_cols <- function(data, overlap_count, keep_gene_seqs = FALSE){
    # get overlapping regions
    names = get_overlap_names(overlap_count)
    data[, paste(names) := get_overlapping_regions(v_gene_sequence, j_gene_sequence, v_trim, j_trim, overlap_count)]

    stopifnot(isFALSE(keep_nt_identity))
    # get MH
    for (pos in c('up', 'down', 'mid')){
        n = paste0('mh_prop_', pos, '_overlap_', overlap_count)
        v_seq_col = paste0('v_gene_', pos, '_overlap_', overlap_count)
        j_seq_col = paste0('j_gene_', pos, '_overlap_', overlap_count)
        cols = c(v_seq_col, j_seq_col)
        subset = unique(data[, ..cols])

        subset[, paste(n) := get_mh_prop(get(v_seq_col), get(j_seq_col))]
        subset[is.na(get(n)), paste(n) := 0]
        data = merge(data, subset, by = cols)
    }

    if (keep_gene_seqs == FALSE) {
        remove = c('v_gene_sequence', 'j_gene_sequence')
        cols = colnames(data)[!(colnames(data) %in% remove)]
        data = data[, ..cols]
    }
    return(data)
}

get_oriented_sequences_for_processed_df <- function(motif_data, gene_type = GENE_NAME){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    for (g in genes){
        temp = copy(whole_nucseq)
        temp_conversion = get_common_genes_from_seqs(temp, g)
        temp_conversion = temp_conversion[substring(get(g), 4, 4) == toupper(substring(g, 1, 1))]
        whole_nucseq = merge(whole_nucseq, temp_conversion, by.x = 'gene', by.y = g, all = TRUE)
    }
    
    cols = c(paste0(genes, '_group'), paste0(genes, '_sequence'))
    gene_cols = c(paste0(genes, '_group'))
    subset_seqs = unique(whole_nucseq[, ..cols])
    subset_seqs[, N := .N, by = gene_cols]

    for (gr in unique(subset_seqs[N > 1]$v_gene_group)){
        temp = subset_seqs[v_gene_group == gr]$v_gene_sequence
        subset_seqs[v_gene_group == gr, v_gene_sequence := temp[1]]
    }
    
    for (gr in unique(subset_seqs[N > 1]$j_gene_group)){
        temp = subset_seqs[j_gene_group == gr]$j_gene_sequence
        subset_seqs[j_gene_group == gr, j_gene_sequence := temp[1]]
    }

    subset_seqs = unique(subset_seqs)

    motif_data = merge(motif_data, subset_seqs[, c('v_gene_group', 'v_gene_sequence')], by = 'v_gene_group')
    motif_data = merge(motif_data, subset_seqs[, c('j_gene_group', 'j_gene_sequence')], by = 'j_gene_group')
    return(motif_data)
}

process_for_mh <- function(motif_data, whole_nucseq = get_oriented_whole_nucseqs(), overlap_vector = c(0, 1, 2, 3, 4), trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    motif_data = get_oriented_sequences_for_processed_df(motif_data)
    
    for (overlap in overlap_vector){
        motif_data = get_mh_prop_cols(motif_data, overlap, keep_gene_seqs = TRUE)
    }

    cols = colnames(motif_data)[!(colnames(motif_data) %like% '_seq')]
    return(motif_data[, ..cols])
}
