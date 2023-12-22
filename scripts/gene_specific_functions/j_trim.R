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

get_all_possible_sites <- function(gene_type = GENE_NAME){
    return(NULL)
}

get_missing_possible_sites <- function(possible_sites, filtered_motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    unobserved = data.table()
    return(unobserved)
}
