get_terminal_sequences <- function(combine_by_terminal){
    # get genes and sequences
    whole_nucseq = get_oriented_whole_nucseqs()
    seqs = whole_nucseq[substring(gene, 4, 4) == toupper(substring(GENE_NAME, 1, 1))]
    
    # group genes by common terminal sequence
    setnames(seqs, 'gene', GENE_NAME)
    map = get_common_genes_from_seqs(seqs)
    together = merge(seqs, map, by = GENE_NAME)

    # get terminal sequences
    terminal_length = UPPER_TRIM_BOUND + REQUIRED_COMMON_NUCS_5 
    together[, terminal_seq := substring(sequence, nchar(sequence)- (terminal_length - 1), nchar(sequence))]
    if (isTRUE(combine_by_terminal)){
        cols = c('terminal_seq', 'gene')
        together = unique(together[, ..cols])
    }
    return(together)
}

get_distances <- function(sequences, combine_by_terminal = TRUE, full_sequence = FALSE, align = FALSE){
    require(ape)
    require(DECIPHER)
    if (isTRUE(full_sequence)){
        seq_list = DNAStringSet(sequences$sequence)
        stopifnot(isFALSE(combine_by_terminal))
        stopifnot(isTRUE(align))
        seq_list = AlignSeqs(seq_list)
    } else {
        seq_list = DNAStringSet(sequences$terminal_seq)
        stopifnot(isFALSE(align))
    }

    if (isTRUE(combine_by_terminal)){
        gene_var = 'gene'
    } else {
        gene_var = GENE_NAME
    }

    dists = DistanceMatrix(seq_list)

    colnames(dists) = sequences[[gene_var]]
    row.names(dists) = sequences[[gene_var]]
    return(dists)
}

get_gene_families <- function(cluster_count, combine_by_terminal = TRUE, full_sequence = FALSE, align = FALSE){
    require(ape)
    require(DECIPHER)
    seqs = get_terminal_sequences(combine_by_terminal)
    dists = get_distances(dists, combine_by_terminal, full_sequence, align)
    dist_format = as.dist(dists) 
    
    clusters = hclust(dist_format)
    clusters_grouped = cutree(clusters, k = cluster_count)
    clusters_grouped_df = as.data.frame(clusters_grouped)
    clusters_grouped_df[[gene_var]] = row.names(clusters_grouped_df)
    
    require(RColorBrewer)
    colors = brewer.pal(cluster_count, 'Set2')

    plot(as.phylo(clusters), type = 'unrooted', tip.color = colors[clusters_grouped], no.margin = TRUE)

    together = merge(seqs, as.data.table(clusters_grouped_df), by = gene_var)
    return(list(cluster_data = together, tree = as.phylo(clusters)))
}

generate_hold_out_sample <- function(motif_data, cluster_genes){
    motif_data_subset = motif_data[!(gene %in% cluster_genes)]
    sample_data = motif_data[gene %in% cluster_genes]
    # updating the total_tcr, p_gene, gene_weight_type, and weighted_observation variables for the newly sampled datasets
    source(paste0('scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'), local = TRUE)
    motif_data_subset = calculate_subject_gene_weight(motif_data_subset)
    source(paste0('scripts/sampling_procedure_functions/p_gene_given_subject.R'), local = TRUE)
    sample_data = calculate_subject_gene_weight(sample_data)
    return(list(sample = sample_data, motif_data_subset = motif_data_subset))
}


