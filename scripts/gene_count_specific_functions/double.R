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
    if ('j_gene_sequence.x' %in% colnames(together)){
        cols = colnames(together)[!(colnames(together) %like% '.x')]
        together = together[, ..cols]
        cols = str_remove(cols, '.y') 
        colnames(together) = cols
    }
    return(together)
}
