get_pwm_score <- function(pwm, motif, positions){
    motif_elements = strsplit(motif, "")[[1]] 
    names(motif_elements) = positions
    score = 0
    for (index in 1:length(motif_elements)){
        element = motif_elements[index]
        element_score = pwm[base == element & parameter == names(element)]$coefficient
        score = score + element_score 
    }
    return(score)
}

get_base_count_pwm_score <- function(pwm, row, trim_type = TRIM_TYPE){
    score = 0
    for (count in c(paste0(trim_type, '_left_base_count_GC'), paste0(trim_type, '_right_base_count_AT'), paste0(trim_type, '_right_base_count_GC'))){
        element = row[[count]]
        element_score = (pwm[parameter == count]$coefficient) * as.numeric(element)
        score = score + element_score
    }
    return(score)
}

get_all_base_count_pwm_score <- function(predicted_trims, pwm, trim_type = TRIM_TYPE){
    counts = c(paste0(trim_type, '_left_base_count_GC'), paste0(trim_type, '_right_base_count_AT'), paste0(trim_type, '_right_base_count_GC'))
    condensed = unique(predicted_trims[,..counts])
    condensed[, pwm_score := apply(condensed, 1, get_base_count_pwm_score, pwm = pwm, trim_type = trim_type)]
    return(condensed)
}

get_all_pwm_score <- function(predicted_trims, pwm, per_gene_resid, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    positions = get_positions(trim_type)
    cols = c(paste0(gene_type, '_group'), trim_type, paste0(trim_type, '_motif'), positions)
    condensed = unique(predicted_trims[,..cols]) 
    # calculate pwm score by gene and trim length
    condensed[, pwm_score := unlist(lapply(get(paste0(trim_type, '_motif')), function(x) as.numeric(get_pwm_score(pwm, x, positions))))]
    if (!is.null(per_gene_resid)){
        tog = merge(condensed, per_gene_resid)
    } else {
        tog = condensed
    }
    return(tog)
}

compare_weight_by_base_count <- function(predicted_trims, trim_type = TRIM_TYPE){
    cols = c(paste0(trim_type, '_left_base_count_GC'), paste0(trim_type, '_right_base_count_AT'), paste0(trim_type, '_right_base_count_GC'))
    avg = predicted_trims[, sum(weighted_observation), by = cols]
    setnames(avg, 'V1', 'total_weight')
    return(avg)
}

compare_weight_by_motif <- function(predicted_trims, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    cols = c(paste0(gene_type, '_group'), trim_type, paste0(trim_type, '_motif'))
    avg = predicted_trims[, sum(weighted_observation), by = cols]
    setnames(avg, 'V1', 'total_weight')
    cols2 = c(paste0(gene_type, '_group'), paste0(trim_type, '_motif'))
    avg2 = avg[, sum(total_weight), by = cols2]
    setnames(avg2, 'V1', 'total_weight_by_gene')
    return(avg2)
}

determine_cluster_count_plot <- function(data, gene_type = GENE_NAME){
    # Determine number of clusters
    data = data[, colSums(data != 0) > 0, with = FALSE]
    scaled = scale(data[, -c(paste0(gene_type, '_group'))])
    wss = (nrow(scaled)-1)*sum(apply(scaled,2,var), na.rm = TRUE)
    for (i in 2:15) wss[i] = sum(kmeans(scaled, centers=i)$withinss)
    plot = plot(1:15, wss, type="b", xlab="Number of Clusters",
    ylab="Within groups sum of squares")
    return(plot)
}

cluster_data <- function(data, cluster_count, cluster_variable, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    data = data[, colSums(data != 0) > 0, with = FALSE]
    scaled = scale(data[, -c(paste0(gene_type, '_group'))])
    # K-Means Cluster Analysis
    fit = kmeans(scaled, cluster_count) 
    # get cluster means 
    aggregate(scaled,by=list(fit$cluster),FUN=mean)
    # append cluster assignment
    scaled = data.table(scaled, fit$cluster)
    scaled[[paste0(gene_type, '_group')]] = data[[paste0(gene_type, '_group')]]
    setnames(scaled, 'V2', 'k_means_cluster')
    # format data 
    long = scaled %>%
        pivot_longer(!c(paste0(gene_type, '_group'), 'k_means_cluster'), names_to = trim_type, values_to = cluster_variable) %>%
        as.data.table()
    return(long)
}

get_pwm_profile_plot_file_path <- function(cluster_count) {
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(MOD_PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0('pwm_profile_', cluster_count, '_clusters')) 
    dir.create(path, recursive = TRUE)
    return(path)
}

predict_trimmming_given_pwm_scores <- function(data, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    data[, exp_pwm_score := exp(pwm_score)]
    col = c(paste0(gene_type, '_group'))
    data[, sum_exp_pwm_score := sum(exp_pwm_score), by = col]
    data[, paste0('predicted_p_', trim_type, '_given_', gene_type) := exp_pwm_score/sum_exp_pwm_score]
    return(data[, -c('exp_pwm_score', 'sum_exp_pwm_score')])
}

get_pwm_prediction_residual_plot_file_path <- function() {
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(MOD_PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, model, paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'PWM_only_prediction_residuals') 
    dir.create(path, recursive = TRUE)
    return(path)
}

get_motif_frequency <- function(motif_data, model_performance_data, motif_of_interest, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    positions = get_positions(trim_type)
    cols = c(paste0(gene_type, '_group'), paste0(trim_type, '_motif'), trim_type, positions)
    motifs = unique(motif_data[, ..cols])
    tog = merge(motifs, model_performance_data)
    col = c(paste0(gene_type, '_group'))
    tog[, freq := .N, by = col]
    cols2 = c(paste0(gene_type, '_group'), paste0(trim_type, '_motif'))
    tog[, paste0(trim_type, '_motif_count') := .N, by = cols2]
    tog = unique(tog[, -c(trim_type)])
    tog[, paste0(trim_type, '_motif_freq') := get(paste0(trim_type, '_motif_count'))/freq]
    tog_subset = tog[get(paste0(trim_type, '_motif')) %in% motif_of_interest]
    other_subset = model_performance_data[!(get(paste0(gene_type, '_group')) %in% tog_subset[[paste0(gene_type, '_group')]])]
    other_subset[, paste0(trim_type, '_motif') := 'CTT/CGT']
    other_subset[, paste0(trim_type, '_motif_freq') := 0]
    tog_subset = rbind(tog_subset, other_subset, fill = TRUE)
    return(tog_subset)
}
