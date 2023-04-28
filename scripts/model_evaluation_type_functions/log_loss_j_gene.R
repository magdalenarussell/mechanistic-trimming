HELD_OUT_FRACTION <<- NA 
REPETITIONS <<- NA 
WRITE_INTERMEDIATE_LOSS <<- NA 

evaluate_loss <- function(motif_data) {
    # fit model to the training data 
    if (MODEL_TYPE != 'null') {
        model = fit_model(motif_data)
        parameter_count = length(model$coefficients) 
    } else {
        model = 'null'
        parameter_count = 0
    } 
    # load J motif data, adjust weighting to reflect actual p_gene_given_subject weight for log loss calculation
    source(paste0(MOD_PROJECT_PATH,'/scripts/sampling_procedure_functions/p_gene_given_subject.R'), local = TRUE)
    v_motif_path = get_subject_motif_output_location() 
    old_GENE_NAME = GENE_NAME
    old_TRIM_TYPE = TRIM_TYPE
    GENE_NAME <<- 'j_gene'
    TRIM_TYPE <<- 'j_trim'
    source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'), local = TRUE)
    j_motif_path = str_replace(v_motif_path, 'v_trim', 'j_trim')
    j_motif_data = aggregate_all_subject_data(j_motif_path)
    j_motif_data = calculate_subject_gene_weight(j_motif_data)
    
    # compute conditional logistic loss value for training data 
    log_loss = calculate_cond_expected_log_loss(model, j_motif_data)

    # reset global variables
    GENE_NAME <<- old_GENE_NAME 
    TRIM_TYPE <<- old_TRIM_TYPE 
    return(list(loss = log_loss, model_parameter_count = parameter_count, held_out_cluster_number = NA, held_out_genes = 'j_genes'))
}
