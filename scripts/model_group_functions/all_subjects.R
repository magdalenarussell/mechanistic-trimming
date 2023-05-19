get_predicted_dist_file_name <- function(subgroup){
    return(paste0(MODEL_GROUP, '_weighting_with_', GENE_WEIGHT_TYPE, '.tsv'))
}

get_pwm_matrix_file_name <- function(subgroup){
    return(paste0(MODEL_GROUP, '_weighting_with_', GENE_WEIGHT_TYPE, '.tsv'))
}

get_coefficient_output_file_name <- function(subgroup){
    return(paste0(MODEL_GROUP, '_weighting_with_', GENE_WEIGHT_TYPE, '.tsv'))
}

fit_model_by_group <- function(motif_data, formula = get_model_formula(trim_type = TRIM_TYPE, gene_type = GENE_NAME), write_coeffs = TRUE, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    model = fit_model(motif_data, formula, trim_type)

    # get predicted probabilities and empirical probabilities
    motif_data$predicted_prob = predict(model, type = 'response')
    cols = c('subject', paste0(gene_type, '_group'))
    motif_data[, empirical_prob := count/sum(count), by = cols]
    motif_data$model_group = MODEL_GROUP

    # save file
    file_path = get_predicted_dist_file_path()
    file_name = get_predicted_dist_file_name(subgroup = NULL)
    location = file.path(file_path, file_name)
   
    if (isTRUE(write_coeffs)){
        fwrite(motif_data, location, sep = '\t')
    }
    
    # get coef path
    coef_path = get_coefficient_output_file_path()
    coef_file_name =  get_coefficient_output_file_name(subgroup = NULL)
    location_coefs = file.path(coef_path, coef_file_name)

    if (grepl('motif', MODEL_TYPE, fixed = TRUE)){
        # calculate coeffiecients
        pwm_matrix = get_coefficient_matrix(motif_data, ref_base = 'A', trim_type = trim_type)$result
        pwm_dt = as.data.table(pwm_matrix)
        pwm_dt$base = rownames(pwm_matrix)
        pwm_dt$model_group = MODEL_GROUP

        # save coefficients
        pwm_file_path = get_pwm_matrix_file_path()
        pwm_file_name = get_pwm_matrix_file_name(subgroup = NULL)
        location = file.path(pwm_file_path, pwm_file_name)
        
        if (isTRUE(write_coeffs)){
            fwrite(pwm_dt, location, sep = '\t')
        }

        # get all coefficients
        coefs = format_model_coefficient_output(model, pwm_dt, trim_type = trim_type)
    } else {
        coefs = format_model_coefficient_output(model, trim_type = trim_type)
    }
    coefs$model_group = MODEL_GROUP
    if (isTRUE(write_coeffs)){
        fwrite(coefs, location_coefs, sep = '\t')
    } else {
        return(coefs)
    }
}

get_predicted_distribution_data <- function(){
    file_path = get_predicted_dist_file_path()
    file_name = get_predicted_dist_file_name(subgroup = NULL)
    location = file.path(file_path, file_name)
 
    predicted_trims = fread(location)
    return(predicted_trims)
}

get_model_coefficient_data <- function(specific_path = NULL){
    if (is.null(specific_path)){
        pwm_file_path = get_coefficient_output_file_path()
    } else {
        pwm_file_path = specific_path
    }
    pwm_file_name = get_pwm_matrix_file_name(subgroup = NULL)
    location = file.path(pwm_file_path, pwm_file_name)
    
    pwm = fread(location)
    return(pwm)
}
