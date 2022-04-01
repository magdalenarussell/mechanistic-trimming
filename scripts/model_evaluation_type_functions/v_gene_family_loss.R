HELD_OUT_FRACTION <<- NA 
REPETITIONS <<- NA 
WRITE_INTERMEDIATE_LOSS <<- NA 

get_terminal_sequences <- function(){
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
    cols = c('terminal_seq', 'gene')
    together = unique(together[, ..cols])
    return(together)
}

get_gene_families <- function(){
    seqs = get_terminal_sequences()
    seq_list = DNAStringSet(seqs$terminal_seq)

    require(DECIPHER)

    dists = DistanceMatrix(seq_list)
    clusts = IdClusters(dists, showPlot = TRUE, cutoff = 0.25)

    colnames(dists) = seqs$gene
    dists = as.data.table(dists)
    dists$gene = seqs$gene

    clusts$gene = seqs$gene

    together = merge(seqs, clusts, by = 'gene')
    together = merge(together, dists, by = 'gene')

    long = together %>%
        pivot_longer(starts_with('TRB'), names_to = 'distance_to_gene', values_to = 'distance') %>%
        as.data.table()

    no_cluster_dists = long
    for (gene_name in unique(long$gene)){
        clust = unique(long[gene == gene_name]$cluster) 
        cluster_genes = unique(long[cluster==clust]$gene)
        no_cluster_dists = no_cluster_dists[!(gene == gene_name & distance_to_gene %in% cluster_genes)] 
    }
    return(no_cluster_dists)
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

evaluate_loss <- function(motif_data) {
    gene_families = get_gene_families()
    gene_families[, min_dist := min(distance), by = cluster]
    held_out_clusters = unique(gene_families[min_dist >= 0.5]$cluster)
    held_out_genes = unique(gene_families[cluster %in% held_out_clusters]$gene)
    log_loss_vector = c()
    parameter_count_vector = c()

    # Generate a held out sample and motif data subset
    sample_data = generate_hold_out_sample(motif_data, held_out_genes) 
    motif_data_subset = sample_data$motif_data_subset
    sample = sample_data$sample

    # Fit model to the motif_data_subset
    if (MODEL_TYPE != 'null'){
        model = fit_model(motif_data_subset)
        parameter_count_vector = c(parameter_count_vector, length(model$coefficients)) 
    } else {
        model = 'null'
        parameter_count_vector = c(parameter_count_vector, 0)
    }

    # Compute conditional logistic loss value for held out sample using model
    log_loss = calculate_cond_expected_log_loss(model, sample)
    log_loss_vector = c(log_loss_vector, log_loss)
       
    held_out_genes = paste(held_out_genes, collapse = ', ')
    held_out_clusters = paste(held_out_clusters, collapse = ', ')
    #TODO: should this be just multiplied by the sample_prob_vector or be the mean? 
    # expected_log_loss = sum(sample_prob_vector * log_loss_vector)

    return(list(loss = log_loss_vector, model_parameter_count = parameter_count_vector, held_out_cluster_number = held_out_clusters, held_out_genes = held_out_genes))
}
