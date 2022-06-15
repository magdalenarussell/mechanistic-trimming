include_snp_in_model_formula <- function(get_model_formula, variable_to_vary_by_snp, snp){
    formula = get_model_formula()
    adapted_formula = str_replace(formula, variable_to_vary_by_snp, paste0('(', variable_to_vary_by_snp, ')*`', snp, '`'))
    adapted_formula[3] = str_replace_all(adapted_formula[3], '\\[', '')
    adapted_formula[3] = str_replace_all(adapted_formula[3], '\\]', '')

    adapted_formula = as.formula(paste0(adapted_formula[2], adapted_formula[1], adapted_formula[3]))
    return(adapted_formula)
}

get_model_snp_analysis_file_path <- function(snp){
    path = file.path(output_path, annotation_type, trim_type, productivity, 'snp_analysis', snp) 
    if (!dir.exists(path)){
        dir.create(path, recursive = true)
    }
    return(path)
}

get_model_snp_analysis_file_name <- function(snp){
    path = get_model_snp_analysis_file_path(snp)
    name = file.path(path, paste0('snp_analysis', model_type, '_', motif_type, '_motif_', left_nuc_motif_count, '_', right_nuc_motif_count, '_bounded_', lower_trim_bound, '_', upper_trim_bound, '_', gene_weight_type, '.tsv'))
    return(name)
}

calculate_lrt_p <- function(original_model_loglik, snp_var_model_loglik, df_diff){
    test_statistic = -2 * (as.numeric(original_model_loglik)-as.numeric(snp_var_model_loglik))
    lrt_p = pchisq(test_statistic, df = df_diff, lower.tail = false)
    return(lrt_p)
}

compile_snp_result <- function(original_model_loglik, snp_var_model_loglik, df_diff, lrt_p, variable_to_vary_by_snp, snp){
    result = data.table(motif_length_5_end = left_nuc_motif_count, motif_length_3_end = right_nuc_motif_count, motif_type = motif_type, gene_weight_type = gene_weight_type, upper_bound = upper_trim_bound, lower_bound = lower_trim_bound, model_type = model_type, terminal_melting_5_end_length = left_side_terminal_melt_length, original_model_loglik = original_model_loglik, variable_to_vary_by_snp = variable_to_vary_by_snp, snp = snp, model_varying_var_by_snp_loglik = snp_var_model_loglik, df_diff = df_diff, likelihood_ratio_p = lrt_p) 
    return(result)
}
 
write_snp_result_dt <- function(result, snp){
    file_name = get_model_snp_analysis_file_name(snp) 
    fwrite(result, file_name, sep = '\t')
    return(result)
}
