get_predicted_dist_file_name <- function(subgroup){
    return(paste0(MODEL_GROUP, '_weighting_with_', GENE_WEIGHT_TYPE, '.tsv'))
}

get_pwm_matrix_file_name <- function(subgroup){
    return(paste0(MODEL_GROUP, '_weighting_with_', GENE_WEIGHT_TYPE, '.tsv'))
}

fit_model_by_group <- function(motif_data, write_coeffs = TRUE){
    model = fit_model(motif_data)

    # get predicted probabilities and empirical probabilities
    motif_data$predicted_prob = predict(model, type = 'response')
    motif_data[, empirical_prob := count/sum(count), by = .(subject, gene)]
    motif_data$model_group = MODEL_GROUP

    # save file
    file_path = get_predicted_dist_file_path()
    file_name = get_predicted_dist_file_name(subgroup = NULL)
    location = file.path(file_path, file_name)
    
    fwrite(motif_data, location, sep = '\t')

    if (grepl('motif', MODEL_TYPE, fixed = TRUE)){
        # calculate coeffiecients
        pwm_matrix = get_coeffiecient_matrix(motif_data, ref_base = 'A')
        pwm_dt = as.data.table(pwm_matrix)
        pwm_dt$base = rownames(pwm_matrix)
        pwm_dt$model_group = MODEL_GROUP

        # save coefficients
        pwm_file_path = get_pwm_matrix_file_path()
        pwm_file_name = get_pwm_matrix_file_name(subgroup = NULL)
        location = file.path(pwm_file_path, pwm_file_name)

        if (isTRUE(write_coeffs)){
            fwrite(pwm_dt, location, sep = '\t')
        } else {
            return(pwm_dt)
        }
    }
}

get_predicted_distribution_data <- function(){
    file_path = get_predicted_dist_file_path()
    file_name = get_predicted_dist_file_name(subgroup = NULL)
    location = file.path(file_path, file_name)
 
    predicted_trims = fread(location)
    return(predicted_trims)
}

get_model_coefficient_data <- function(){
    pwm_file_path = get_pwm_matrix_file_path()
    pwm_file_name = get_pwm_matrix_file_name(subgroup = NULL)
    location = file.path(pwm_file_path, pwm_file_name)
    
    pwm = fread(location)
    return(pwm)
}
