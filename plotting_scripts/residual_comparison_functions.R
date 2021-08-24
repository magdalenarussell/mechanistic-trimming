if (!is.null(RESIDUAL_COMPARE_FEATURE)){
    source(paste0('plotting_scripts/residual_comparison_features/', RESIDUAL_COMPARE_FEATURE, '.R'))
}

# calculate root mean square error
calculate_rmse <- function(predicted_trims){
    total_trims = unique(predicted_trims$trim_length)
    # calculate subjects per gene
    subj_counts = predicted_trims[, .N/length(total_trims), by = gene]
    setnames(subj_counts, 'V1', 'subject_count')
    
    # calculate residuals for each subject, gene, trim
    predicted_trims[, residual := empirical_prob - predicted_prob]
    predicted_trims[, residual_sq := (residual)^2]

    # calculate residual sum by subject, gene
    resid_sq_sums = predicted_trims[, sum(residual_sq), by = .(gene, subject)]
    setnames(resid_sq_sums, 'V1', 'resid_sq_sum')

    # calculate mean residual across trims
    together = merge(resid_sq_sums, subj_counts)
    together[, mean_resid_sq := resid_sq_sum/length(total_trims)]

    # calculate mean residual sum by gene, calculate root mean square error by dividing by subject count
    mean_resid_sq_sums = together[, sum(mean_resid_sq), by = .(gene, subject_count)]
    setnames(mean_resid_sq_sums, 'V1', 'mean_resid_sq_sum')
    mean_resid_sq_sums[, rmse := sqrt(mean_resid_sq_sum/subject_count)]

    return(mean_resid_sq_sums[, -c('mean_resid_sq_sum')])
}

# calculate root mean square error
calculate_residual_by_position <- function(predicted_trims){
    # calculate residuals for each subject, gene, trim
    predicted_trims[, residual := empirical_prob - predicted_prob]

    avg_resids = predicted_trims[, mean(residual), by = .(gene, trim_length)]
    setnames(avg_resids, 'V1', 'avg_resid')
    return(avg_resids)
}
