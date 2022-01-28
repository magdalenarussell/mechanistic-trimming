HELD_OUT_FRACTION <<- NA 
REPETITIONS <<- NA 
WRITE_INTERMEDIATE_LOSS <<- NA 

generate_hold_out_sample_old_loss_cv <- function(motif_data, sample_size){
     stopifnot(sample_size < length(unique(motif_data$subject))*length(unique(motif_data$gene)))
     # sample size is the number of gene, subject combos
     set.seed(66)
     sample_genes = sample(unique(motif_data$gene), sample_size, replace = TRUE)
     sample_data = data.table()
     motif_data_subset = motif_data
     for (sample in sample_genes){
        sample_subject = sample(unique(motif_data$subject), 1)
        temp = motif_data[subject == sample_subject & gene == sample]
        sample_data = rbind(sample_data,temp)
        cols = colnames(sample_data)
        motif_data_subset = motif_data_subset[!sample_data, on = cols]
     }
     return(list(sample = sample_data, motif_data_subset = motif_data_subset))
}

evaluate_loss <- function(motif_data) {
    # Generate a held out sample and motif data subset
    sample_data = generate_hold_out_sample_old_loss_cv(motif_data, sample_size = length(unique(motif_data$gene))) 
    motif_data_subset = sample_data$motif_data_subset
    sample = sample_data$sample

    # Fit model to the motif_data_subset
    model = fit_model(motif_data_subset)
    parameter_count = length(model$coefficients)

    # Compute conditional logistic loss value for held out sample using model
    source(paste0('scripts/sampling_procedure_functions/raw_count.R'), local = TRUE)
    sample = calculate_subject_gene_weight(sample)
    log_loss = calculate_cond_expected_log_loss(model, sample)

    return(list(loss = log_loss, model_parameter_count = unique(parameter_count)))
}
