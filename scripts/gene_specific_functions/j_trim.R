JOINING_GENE <<- 'd_gene'
JOINING_TRIM <<- 'd1_trim'

get_oriented_whole_nucseqs <- function(whole_nucseq = get_whole_nucseqs()){
    # reorient sequence so that it is 5 -> 3 on the actual trimmed strand
    whole_nucseq[, sequence := unlist(lapply(sequence, function(x) as.character(reverseComplement(DNAString(x)))))]
    return(whole_nucseq)
}

