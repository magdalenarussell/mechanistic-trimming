get_oriented_whole_nucseqs <- function(whole_nucseq = get_whole_nucseqs()){
    # reorient sequence so that it is 5 -> 3 on the actual trimmed strand
    whole_nucseq[, sequence := unlist(lapply(sequence, function(x) as.character(reverseComplement(DNAString(x)))))]
    setnames(whole_nucseq, 'sequence', paste0('j_gene_sequence'))
    return(whole_nucseq)
}

