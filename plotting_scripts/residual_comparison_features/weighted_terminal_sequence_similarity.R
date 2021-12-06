
get_feature <- function(predicted_trims){
    whole_nucseq = fread(get(paste0('WHOLE_NUCSEQS_', ANNOTATION_TYPE)))
    gene_seqs = whole_nucseq[substring(gene_names, 4, 4) == toupper(substring(GENE_NAME, 1, 1))]
    colnames(gene_seqs) = c(GENE_NAME, 'sequences')
    gene_groups = get_common_genes_from_seqs(gene_seqs)
    together = merge(gene_groups, gene_seqs, by = 'v_gene')
    together$terminal_seqs = get_terminal_seq(together$sequences)

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
    motif_length = RIGHT_NUC_MOTIF_COUNT + LEFT_NUC_MOTIF_COUNT
    weights_left = seq(1, motif_length)
    weights_right = seq(motif_length, 1)
    remaining_weights = rep(motif_length, nchar(ref_seq) - length(weights_left) - length(weights_right))
    weights = c(weights_left, remaining_weights, weights_right)
    stopifnot(length(weights) == nchar(ref_seq))

    for (seq in terminal_seq_list){
        dist = 0
        for (base in seq(nchar(ref_seq))){
            if (substring(ref_seq, base, base) != substring(seq, base, base)){
                position_value = weights[base]
                dist = dist + position_value 
            } 
        }
        dists = c(dists, dist)
    }
    return(dists)
}

