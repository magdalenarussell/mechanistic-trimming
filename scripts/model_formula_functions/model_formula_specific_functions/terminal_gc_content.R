get_GC_content <- function(){
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

    # get length of terminal seq
    together[, depth := trim_length + LEFT_NUC_MOTIF_COUNT]

    # get terminal seq (only interested in double stranded end (not including pnucs or pnuc pairs))
    together[, terminal_seq := substring(sequence, nchar(sequence) - depth + 1, nchar(sequence)-abs(PNUC_COUNT))]

    # get GC content
    seq_list = DNAStringSet(together$terminal_seq)
    base_counts = as.data.table(letterFrequency(seq_list, letters="ACGT", OR = 0))
    base_counts[, terminal_gc_content := (C+G)/(A+C+G+T)]

    # merge
    together = cbind(together, base_counts)
    return(unique(together[, c('gene', 'trim_length', 'terminal_gc_content')]))
}

process_for_terminal_gc_content <- function(group_motif_data){
    row_count = nrow(group_motif_data)
    if (!('terminal_gc_content' %in% colnames(group_motif_data))){
        terminal_gc = get_GC_content()
        if (is.factor(group_motif_data$trim_length)){ 
            terminal_gc$trim_length = as.factor(terminal_gc$trim_length)
        }
        together = merge(group_motif_data, terminal_gc, by = c('gene', 'trim_length'))
    } else {
        together = group_motif_data
    }
    stopifnot(nrow(together) == row_count)
    return(together)
}

