HELD_OUT_FRACTION <<- NA 
REPETITIONS <<- NA 
WRITE_INTERMEDIATE_LOSS <<- NA 

evaluate_loss <- function(motif_data, trim_type = TRIM_TYPE, gene_type = GENE_NAME) {
    # fit model to the training data 
    if (MODEL_TYPE != 'null') {
        model = fit_model(motif_data, trim_type= trim_type)
        parameter_count = length(model$coefficients) 
    } else {
        model = 'null'
        parameter_count = 0
    } 
    # load J motif data, adjust weighting to reflect actual p_gene_given_subject weight for log loss calculation
    source(paste0(MOD_PROJECT_PATH,'/scripts/sampling_procedure_functions/p_gene_pooled.R'), local = TRUE)
    v_motif_path = get_subject_motif_output_location() 
    old_GENE_NAME = GENE_NAME
    old_TRIM_TYPE = TRIM_TYPE
    GENE_NAME <<-  JOINING_GENE
    TRIM_TYPE <<- JOINING_TRIM
    source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'), local = TRUE)
    j_motif_path = str_replace(v_motif_path, old_TRIM_TYPE, TRIM_TYPE)
    j_motif_data = aggregate_all_subject_data(j_motif_path, trim_type = TRIM_TYPE)
    j_motif_data = calculate_subject_gene_weight(j_motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE)
    
    # compute conditional logistic loss value for training data 
    log_loss = calculate_cond_expected_log_loss(model, j_motif_data, trim_type = trim_type, gene_type = gene_type, joining_trim = JOINING_TRIM, joining_gene_type = JOINING_GENE)

    # reset global variables
    GENE_NAME <<- old_GENE_NAME 
    TRIM_TYPE <<- old_TRIM_TYPE 
    return(list(loss = log_loss, model_parameter_count = parameter_count, held_out_cluster_number = NA, held_out_genes = JOINING_GENE))
}
