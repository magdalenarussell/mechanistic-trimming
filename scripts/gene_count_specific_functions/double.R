get_oriented_full_sequences <- function(subject_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME){
    genes = get_gene_order(gene_type)
    gene1seq = paste0(genes[1], '_sequence')
    gene2seq = paste0(genes[2], '_sequence')

    temp_data = merge(subject_data, whole_nucseq[,!..gene2seq], by.x = genes[1], by.y = 'gene')
    temp_data = merge(temp_data, whole_nucseq[,!..gene1seq], by.x = genes[2], by.y = 'gene')
    together = temp_data
    for (g in genes) {
        gene_seqs = whole_nucseq[substring(gene, 4, 4) %in% toupper(substring(g, 1, 1))]
        gene_groups = get_common_genes_from_seqs(gene_seqs, g)
        together = merge(together, gene_groups, by = g)
    }
    return(together)
}

# TODO complete this...this should only be necessary if INSERTIONS == FALSE
filter_motif_data_for_possible_sites <- function(motif_data){
    # if (JOINING_INSERT == 'vj_insert'){
        # frame_data = get_frames_data()
    # }
    return(motif_data)
}
