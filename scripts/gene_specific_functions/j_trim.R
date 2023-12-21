get_oriented_whole_nucseqs <- function(whole_nucseq = get_whole_nucseqs()){
    # reorient sequence so that it is 5 -> 3 on the actual trimmed strand
    whole_nucseq[, sequence := unlist(lapply(sequence, function(x) as.character(reverseComplement(DNAString(x)))))]
    setnames(whole_nucseq, 'sequence', paste0('j_gene_sequence'))
    return(whole_nucseq)
}

source(paste0(MOD_PROJECT_PATH,'/scripts/gene_count_specific_functions/single.R'))

filter_motif_data_for_possible_sites <- function(motif_data){
    return(motif_data)
}

adjust_trimming_sites_for_ligation_mh <- function(motif_data){
    return(motif_data)
}
