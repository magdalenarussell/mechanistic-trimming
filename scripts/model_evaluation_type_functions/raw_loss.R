HELD_OUT_FRACTION <<- NA 
REPETITIONS <<- NA 
WRITE_INTERMEDIATE_LOSS <<- NA 

evaluate_loss <- function(motif_data){
    # Fit model to the training data 
    model = fit_model(motif_data)
    parameter_count = length(model$coefficients) 
    # Adjust weighting to reflect actual raw counts for log loss calculation
    source(paste0('scripts/sampling_procedure_functions/raw_count.R'), local = TRUE)
    motif_data = calculate_subject_gene_weight(motif_data)
    # Compute conditional logistic loss value for training data 
    log_loss = calculate_cond_expected_log_loss(model, motif_data)
    return(list(loss = log_loss, model_parameter_count = parameter_count, held_out_cluster_number = NA, held_out_genes = NA))
}
