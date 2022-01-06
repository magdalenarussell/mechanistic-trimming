get_pwm_score <- function(pwm, motif, positions){
    motif_elements = strsplit(motif, "")[[1]] 
    names(motif_elements) = positions
    score = 0
    for (index in 1:length(motif_elements)){
        element = motif_elements[index]
        element_score = pwm[base == element][[names(element)]]
        score = score + element_score 
    }
    return(score)
}

determine_cluster_count_plot <- function(data){
    # Determine number of clusters
    data = data[, colSums(data != 0) > 0, with = FALSE]
    scaled = scale(data[, -c('gene')])
    wss = (nrow(scaled)-1)*sum(apply(scaled,2,var), na.rm = TRUE)
    for (i in 2:15) wss[i] = sum(kmeans(scaled, centers=i)$withinss)
    plot = plot(1:15, wss, type="b", xlab="Number of Clusters",
    ylab="Within groups sum of squares")
    return(plot)
}

cluster_data <- function(data, cluster_count, cluster_variable){
    data = data[, colSums(data != 0) > 0, with = FALSE]
    scaled = scale(data[, -c('gene')])
    # K-Means Cluster Analysis
    fit = kmeans(scaled, cluster_count) 
    # get cluster means 
    aggregate(scaled,by=list(fit$cluster),FUN=mean)
    # append cluster assignment
    scaled = data.table(scaled, fit$cluster)
    scaled$gene = data$gene
    setnames(scaled, 'V2', 'k_means_cluster')
    # format data 
    long = scaled %>%
        pivot_longer(!c('gene', 'k_means_cluster'), names_to = 'trim_length', values_to = cluster_variable) %>%
        as.data.table()
    return(long)
}

get_pwm_profile_plot_file_path <- function(cluster_count) {
    if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0('pwm_profile_', cluster_count, '_clusters')) 
    dir.create(path, recursive = TRUE)
    return(path)
}

predict_trimmming_given_pwm_scores <- function(data){
    data[, exp_pwm_score := exp(pwm_score)]
    data[, sum_exp_pwm_score := sum(exp_pwm_score), by = gene]
    data[, predicted_p_trim_given_gene := exp_pwm_score/sum_exp_pwm_score]
    return(data[, -c('exp_pwm_score', 'sum_exp_pwm_score')])
}

get_pwm_prediction_residual_plot_file_path <- function() {
    if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'PWM_only_prediction_residuals') 
    dir.create(path, recursive = TRUE)
    return(path)
}


