get_feature <- function(predicted_trims){
    whole_nucseq = get_oriented_whole_nucseqs()
    gene_seqs = whole_nucseq[substring(gene, 4, 4) == toupper(substring(GENE_NAME, 1, 1))]
    setnames(gene_seqs, 'gene', GENE_NAME)
    colnames(gene_seqs) = c(GENE_NAME, 'sequence')
    gene_groups = get_common_genes_from_seqs(gene_seqs)
    together = merge(gene_groups, gene_seqs, by = GENE_NAME)
    together$terminal_seqs = get_terminal_seq(together$sequence)
    together$GC_freq = get_GC_content(together$terminal_seqs) 
    subset = together[, c('gene', 'GC_freq')]
    setnames(subset, 'GC_freq', 'feature')
    return(subset)
}

get_terminal_seq <- function(whole_seq_list){
    terminal_seqs = c()
    for (index in seq(1, length(whole_seq_list))){
        temp_seq = DNAStringSet(whole_seq_list[index])
        length = nchar(temp_seq)
        term_length = UPPER_TRIM_BOUND + LEFT_NUC_MOTIF_COUNT
        term_seq = substring(temp_seq, length - term_length + 1, length)
        terminal_seqs = c(terminal_seqs, term_seq)
    }
    return(terminal_seqs)
}

get_GC_content <- function(seq_list){
    seq_list = DNAStringSet(seq_list)
    base_counts = as.data.table(letterFrequency(seq_list, letters="ACGT", OR = 0))
    base_counts[, GC_freq := (C+G)/(A+C+G+T)]
    return(base_counts$GC_freq)
}

