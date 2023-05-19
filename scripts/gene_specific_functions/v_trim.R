get_oriented_whole_nucseqs <- function(whole_nucseq = get_whole_nucseqs()){
    setnames(whole_nucseq, 'sequence', 'v_gene_sequence')
    return(whole_nucseq)
}

