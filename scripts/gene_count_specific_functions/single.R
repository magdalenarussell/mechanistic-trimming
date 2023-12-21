get_oriented_full_sequences <- function(subject_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME){
    temp_data = merge(subject_data, whole_nucseq, by.x = gene_type, by.y = 'gene')
    gene_seqs = whole_nucseq[substring(gene, 4, 4) == toupper(substring(gene_type, 1, 1))]
    gene_groups = get_common_genes_from_seqs(gene_seqs, gene_type)
    together = merge(temp_data, gene_groups, by = gene_type)
    return(together)
}
