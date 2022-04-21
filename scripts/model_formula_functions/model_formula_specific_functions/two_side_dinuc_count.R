get_all_dinucs_variables <- function(side){
    stopifnot(side %in% c('left', 'right', 'span'))
    if (side %in% c('left', 'right')) {
        bases = c('AT_GC', 'GC_AT', 'GC_GC', 'AT_AT')
        vars = paste0(side, '_dinuc_count_', bases)
    } else {
        vars = paste0(side, '_identity')
    }
    return(vars)
}

convert_dinucs_to_identity_dinucs <- function(dinuc_vector){
    nuc_identity = list('A' = 'AT', 'T' = 'AT', 'G' = 'GC', 'C' = 'GC')
    identity_dinucs = c() 
    for (dinuc in dinuc_vector){
        first_base = substring(dinuc, 1, 1)
        second_base = substring(dinuc, 2, 2)
        id_dinuc = paste0(nuc_identity[[first_base]], '_', nuc_identity[[second_base]])
        identity_dinucs = c(identity_dinucs, id_dinuc)
    }
    together = data.table(span_seq = dinuc_vector, span_identity = identity_dinucs)
    return(together)
}

count_dinucs_seq_list <- function(seq_list, side){
    subset = unique(seq_list)
    bases = c('AT_AT', 'AT_GC', 'GC_AT', 'GC_GC')
    registerDoParallel(cores=NCPU)  
    cl = makeCluster(NCPU, type="FORK")
    seq_list_DNA = parLapply(cl, subset, function(x) DNAString(x))
    counts = parLapply(cl, seq_list_DNA, function(x) dinucleotideFrequency(x))
    counts_dt = rbindlist(lapply(counts, as.data.frame.list))
    counts_dt[, 'AT_AT' := AA+AT+TA+TT]
    counts_dt[, 'GC_AT' := GA+GT+CA+CT]
    counts_dt[, 'AT_GC' := AG+AC+TG+TC]
    counts_dt[, 'GC_GC' := GG+GC+CG+CC]
    stopCluster(cl)  
    paired_cols = c()
    for (base in bases) {
        new_name = paste0(side, '_dinuc_count_', base)
        setnames(counts_dt, base, new_name)
        paired_cols = c(paired_cols, new_name)
    }
    together = cbind(subset, counts_dt[, ..paired_cols])
    setnames(together, 'subset', paste0(side, '_seq'))
    return(together)
}
    
get_left_right_seq_vars_dinuc <- function(motif_data, left_nuc_count = LEFT_SIDE_TERMINAL_MELT_LENGTH){
    whole_nucseq = get_oriented_whole_nucseqs()
    trims = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)
    
    genes = whole_nucseq$gene[substring(whole_nucseq$gene, 4, 4) == toupper(substring(GENE_NAME, 1, 1))]
    together = data.table(gene = rep(genes, length(trims)), trim_length = rep(trims, length(genes)))
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
    together[, right_seq := substring(sequence, nchar(sequence) - trim_length + 1, nchar(sequence)-abs(PNUC_COUNT))]
    together[, span_seq := substring(sequence, nchar(sequence) - trim_length, nchar(sequence)- trim_length + 1)]

    if (is.numeric(left_nuc_count)){
        together[, left_seq := substring(sequence, nchar(sequence) - (trim_length + left_nuc_count) + 1, nchar(sequence)-trim_length)]
    } else if (left_nuc_count == 'right_nuc_count') {
        together[, left_seq := substring(sequence, nchar(sequence) - (trim_length + nchar(right_seq)) + 1, nchar(sequence)-trim_length)]
    }

    motif_data_together = merge(motif_data, unique(together[, -c('sequence')]), by = c('gene', 'trim_length'))
    return(motif_data_together)
}

process_for_two_side_dinuc_count <- function(motif_data, left_nuc_count = LEFT_SIDE_TERMINAL_MELT_LENGTH){
    vars = get_all_dinucs_variables('left')
    if (!all(vars %in% colnames(motif_data))){
        motif_data = get_left_right_seq_vars_dinuc(motif_data, left_nuc_count)

        for (side in c('left', 'right')){
            col = paste0(side, '_seq')
            dinuc_counts = count_dinucs_seq_list(motif_data[[col]], side)
            motif_data = merge(motif_data, dinuc_counts, by = col)
        }
        span = convert_dinucs_to_identity_dinucs(unique(motif_data$span_seq))
        motif_data = merge(motif_data, span, by = 'span_seq')
    }
    return(motif_data[, -c('left_seq', 'right_seq')])
}
