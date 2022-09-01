get_all_base_variables <- function(side){
    stopifnot(side %in% c('left', 'right'))
    bases = c('GC', 'AT')
    vars = paste0(side, '_base_count_', bases)
    return(vars)
}

count_bases_seq_list <- function(seq_list, side){
    subset = unique(seq_list)
    bases = c('A', 'T', 'G', 'C')
    registerDoParallel(cores=NCPU)  
    cl = makeCluster(NCPU, type="FORK")
    seq_list_DNA = parLapply(cl, subset, function(x) DNAString(x))
    counts = parLapply(cl, seq_list_DNA, function(x) letterFrequency(x, letters = bases))
    counts_dt = rbindlist(lapply(counts, as.data.frame.list))
    counts_dt[, 'AT' := A+T]
    counts_dt[, 'GC' := G+C]
    stopCluster(cl)  
    paired_cols = c()
    for (base in c('AT', 'GC')) {
        new_name = paste0(side, '_base_count_', base)
        setnames(counts_dt, base, new_name)
        paired_cols = c(paired_cols, new_name)
    }
    together = cbind(subset, counts_dt[, ..paired_cols])
    setnames(together, 'subset', paste0(side, '_seq'))
    return(together)
}
    
get_left_right_seq_vars <- function(motif_data, left_nuc_count = LEFT_SIDE_TERMINAL_MELT_LENGTH, beyond_motif, single_stranded = FALSE, whole_nucseq = get_oriented_whole_nucseqs()){
    trims = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)
    
    genes = whole_nucseq$gene[substring(whole_nucseq$gene, 4, 4) == toupper(substring(GENE_NAME, 1, 1))]
    together = data.table(gene = rep(genes, each = length(trims)), trim_length = rep(trims, length(genes)))
    together = merge(together, whole_nucseq, by = 'gene')

    setnames(together, 'gene', GENE_NAME)
    map = get_common_genes_from_seqs(together)
    together = merge(together, map, by = GENE_NAME)
    cols = c('trim_length', 'sequence', 'gene')
    together = together[, ..cols]

    if ('name' %in% colnames(together)){
        together = together[, -c('name')]
    }

    # get terminal seq
    if (isTRUE(beyond_motif)){
        right_position_shift = RIGHT_NUC_MOTIF_COUNT
        left_position_shift = LEFT_NUC_MOTIF_COUNT 
    } else {
        right_position_shift = 0
        left_position_shift = 0
    }

    if (isFALSE(single_stranded)){
        together[, right_seq := substring(sequence, nchar(sequence) - trim_length + 1 + right_position_shift, nchar(sequence)-abs(PNUC_COUNT))]
    } else {
        together[, right_seq := substring(sequence, nchar(sequence) - trim_length + 1 + right_position_shift, nchar(sequence))]
    }

    if (is.numeric(left_nuc_count)){
        together[, left_seq := substring(sequence, nchar(sequence) - (trim_length + left_nuc_count) + 1, nchar(sequence)-trim_length - left_position_shift)]
    } else if (left_nuc_count == 'right_nuc_count') {
        together[, left_seq := substring(sequence, nchar(sequence) - (trim_length + nchar(right_seq)) + 1, nchar(sequence)-trim_length - left_position_shift)]
    }
   
    motif_data_together = merge(motif_data, unique(together[, -c('sequence')]), by = c('gene', 'trim_length'))
    return(motif_data_together)
}

process_for_two_side_base_count <- function(motif_data, left_nuc_count = LEFT_SIDE_TERMINAL_MELT_LENGTH, beyond_motif = FALSE, single_stranded = FALSE, whole_nucseq = get_oriented_whole_nucseqs()){
    vars = get_all_base_variables('left')
    if (!all(vars %in% colnames(motif_data))){
        motif_data = get_left_right_seq_vars(motif_data, left_nuc_count, beyond_motif, single_stranded, whole_nucseq)

        for (side in c('left', 'right')){
            col = paste0(side, '_seq')
            base_counts = count_bases_seq_list(motif_data[[col]], side)
            motif_data = merge(motif_data, base_counts, by = col)
        }
    }
    return(motif_data[, -c('left_seq', 'right_seq')])
}
