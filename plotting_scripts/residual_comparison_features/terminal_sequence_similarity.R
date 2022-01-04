get_feature <- function(predicted_trims){
    whole_nucseq = get_oriented_whole_nucseqs()
    gene_seqs = whole_nucseq[substring(gene, 4, 4) == toupper(substring(GENE_NAME, 1, 1))]
    setnames(gene_seqs, 'gene', GENE_NAME)
    colnames(gene_seqs) = c(GENE_NAME, 'sequence')
    gene_groups = get_common_genes_from_seqs(gene_seqs)
    together = merge(gene_groups, gene_seqs, by = GENE_NAME)
    together$terminal_seqs = get_terminal_seq(together$sequence)

    resids = calculate_rmse(predicted_trims)
    smallest_resid = resids[order(rmse)][1] 
    ref_seq = together[gene == smallest_resid$gene]$terminal_seq

    together$feature = get_seq_similarity(ref_seq, together$terminal_seq)
    subset = together[, c('gene', 'feature')]
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

get_seq_similarity <- function(ref_seq, terminal_seq_list){
    dists = c()
    for (seq in terminal_seq_list){
        dist = stringDist(c(ref_seq, seq), method = 'hamming')
        dists = c(dists, dist)
    }
    return(dists)
}

