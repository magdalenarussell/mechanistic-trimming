JOINING_GENE <<- 'd_gene'
JOINING_TRIM <<- 'd0_trim'

if (ANNOTATION_TYPE %like% 'alpha'){
    JOINING_GENE <<- 'j_gene'
    JOINING_TRIM <<- 'j_trim'
}

get_oriented_whole_nucseqs <- function(whole_nucseq = get_whole_nucseqs()){
    return(whole_nucseq)
}

