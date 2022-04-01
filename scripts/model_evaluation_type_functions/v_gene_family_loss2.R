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

get_gene_families <- function(cluster_count){
    seqs = get_terminal_sequences()
    seq_list = DNAStringSet(seqs$terminal_seq)

    require(ape)
    require(DECIPHER)

    dists = DistanceMatrix(seq_list)

    colnames(dists) = seqs$gene
    row.names(dists) = seqs$gene
    dist_format = as.dist(dists) 
    
    clusters = hclust(dist_format)
    clusters_grouped = cutree(clusters, k = cluster_count)
    clusters_grouped_df = as.data.frame(clusters_grouped)
    clusters_grouped_df$gene = row.names(clusters_grouped_df)
    
    require(RColorBrewer)
    colors = brewer.pal(cluster_count, 'Set2')

    plot(as.phylo(clusters), type = 'unrooted', tip.color = colors[clusters_grouped], no.margin = TRUE)

    together = merge(seqs, as.data.table(clusters_grouped_df), by = 'gene')
    return(together)
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
    gene_families = get_gene_families(cluster_count = 3)
    cluster_counts = gene_families[, .N, by = clusters_grouped]
    largest_cluster = cluster_counts[N == max(N)]$clusters_grouped
    held_out_clusters = unique(gene_families[clusters_grouped != largest_cluster]$clusters_grouped)
    log_loss_vector = c()
    parameter_count_vector = c()
    genes = c()
    clusters = c()
    for (cluster_group in c(as.list(held_out_clusters), list(held_out_clusters))) {
        held_out_genes = unique(gene_families[clusters_grouped %in% cluster_group]$gene)
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
        genes = c(genes, held_out_genes)
        cluster_string = paste(cluster_group, collapse = ', ')
        clusters = c(clusters, cluster_string)
    }
    return(list(loss = log_loss_vector, model_parameter_count = parameter_count_vector, held_out_cluster_number = clusters, held_out_genes = genes))
}
