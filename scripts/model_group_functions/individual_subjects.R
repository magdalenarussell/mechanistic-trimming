get_predicted_dist_file_name <- function(subgroup){
    return(paste0(subgroup, '_', MODEL_GROUP, '_weighting_with_', GENE_WEIGHT_TYPE, '.tsv'))
}

get_pwm_matrix_file_name <- function(subgroup){
    return(paste0(subgroup, '_', MODEL_GROUP, '_weighting_with_', GENE_WEIGHT_TYPE, '.tsv'))
}

get_coefficient_output_file_name <- function(subgroup){
    return(paste0(subgroup, '_', MODEL_GROUP, '_weighting_with_', GENE_WEIGHT_TYPE, '.tsv'))
}


fit_model_by_group <- function(motif_data){
    options(contrasts = rep("contr.sum", 2))

    for (indiv in unique(motif_data$subject)){
        motif_data_subset = motif_data[subject == indiv]
        model = fit_model(motif_data_subset)

        # get predicted probabilities and empirical probabilities
        motif_data_subset$predicted_prob = predict(model, type = 'response')
        motif_data_subset[, empirical_prob := count/sum(count), by = .(subject, gene)]
        motif_data_subset$model_group = indiv 

        # save file
        file_path = get_predicted_dist_file_path()
        file_name = get_predicted_dist_file_name(subgroup = indiv)
        complete_name = file.path(file_path, file_name)
        fwrite(motif_data_subset, complete_name, sep = '\t')

        # get coef path
        coef_path = get_coefficient_output_file_path()
        coef_file_name =  get_coefficient_output_file_name(subgroup = indiv)
        location_coefs = file.path(coef_path, coef_file_name)

        if (grepl('motif', MODEL_TYPE, fixed = TRUE)){
            # calculate coeffiecients
            pwm_matrix = get_coeffiecient_matrix(motif_data_subset, ref_base = 'A')
            pwm_dt = as.data.table(pwm_matrix)
            pwm_dt$base = rownames(pwm_matrix)
            pwm_dt$model_group = indiv

            # save coefficients
            pwm_file_path = get_pwm_matrix_file_path()
            pwm_file_name = get_pwm_matrix_file_name(subgroup = indiv)
            pwm_complete_name = file.path(pwm_file_path, pwm_file_name)
            fwrite(pwm_dt, pwm_complete_name, sep = '\t')
            # get all coefficients
            coefs = format_model_coefficient_output(model, pwm_dt)
        } else {
            coefs = format_model_coefficient_output(model)
        }
        coefs$model_group = indiv
        fwrite(coefs, location_coefs, sep = '\t')
    }
}

get_predicted_distribution_data <- function(){
    file_path = get_predicted_dist_file_path()
    files = fs::dir_ls(path = file_path)
    registerDoParallel(cores=NCPU)
    together = foreach(file = files, .combine=rbind) %dopar% {
        fread(file)
    }
    stopImplicitCluster()

    return(together)
}

get_model_coefficient_data <- function(){
    pwm_file_path = get_pwm_matrix_file_path()
    files = fs::dir_ls(path = pwm_file_path)
    registerDoParallel(cores=NCPU)
    together = foreach(file = files, .combine=rbind) %dopar% {
        fread(file)
    }
    stopImplicitCluster()

    return(together)
}

get_all_subject_model_coefficient_data <- function(){
    pwm_file_path = get_pwm_matrix_file_path()
    pwm_file_name = get_pwm_matrix_file_name(subgroup = NULL)
    pwm_file_name = str_sub(pwm_file_name, 2)
    location = file.path(pwm_file_path, pwm_file_name)
    location = str_replace_all(location, MODEL_GROUP, 'all_subjects') 

    pwm = fread(location)
    return(pwm)
}
