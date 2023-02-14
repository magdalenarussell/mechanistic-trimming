HELD_OUT_FRACTION <<- 0.4
REPETITIONS <<- 20

generate_hold_out_sample <- function(motif_data, sample_size){
    stopifnot(sample_size < length(unique(motif_data$gene)))
    # sample size is the number of gene, subject combos
    sample_genes = sample(unique(motif_data$gene), sample_size, replace = FALSE)
    sample_data = data.table()
    motif_data_subset = motif_data
    for (sample in sample_genes){
        temp = motif_data[gene == sample]
        sample_data = rbind(sample_data, temp)
        cols = colnames(sample_data)
        motif_data_subset = motif_data_subset[!sample_data, on = cols]
    }
    # updating the total_tcr, p_gene, gene_weight_type, and weighted_observation variables for the newly sampled datasets
    source(paste0('scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'), local = TRUE)
    motif_data_subset = calculate_subject_gene_weight(motif_data_subset)
    source(paste0('scripts/sampling_procedure_functions/p_gene_given_subject.R'), local = TRUE)
    sample_data = calculate_subject_gene_weight(sample_data)
    return(list(sample = sample_data, motif_data_subset = motif_data_subset))
}

get_hold_out_sample_probability <- function(sample_size, motif_data){
    total_genes = length(unique(motif_data$gene))
    prob = (1/total_genes)^sample_size
    return(prob)
}

evaluate_loss <- function(motif_data) {
    set.seed(66)
    gene_count = length(unique(motif_data$gene))
    held_out_gene_count = round(HELD_OUT_FRACTION*gene_count)
    log_loss_vector = c()
    sample_prob_vector = c()
    parameter_count_vector = c()
    for (rep in 1:REPETITIONS){
        # Generate a held out sample and motif data subset
        sample_data = generate_hold_out_sample(motif_data, sample_size = held_out_gene_count) 
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

        # Compute probability of held out sample
        prob = get_hold_out_sample_probability(held_out_gene_count, motif_data)
        sample_prob_vector = c(sample_prob_vector, prob)
    }

    expected_log_loss = sum((1/REPETITIONS) * log_loss_vector)

    return(list(loss = expected_log_loss, model_parameter_count = unique(parameter_count_vector), held_out_cluster_number = NA, held_out_genes = 'averaged', vect = log_loss_vector))
}
