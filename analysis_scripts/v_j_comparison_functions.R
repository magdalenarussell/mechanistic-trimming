include_gene_type_in_model_formula <- function(get_model_formula, variable_to_vary){
    formula = get_model_formula()
    adapted_formula = str_replace(formula, variable_to_vary, paste0('(', variable_to_vary, ')* gene_type'))
    adapted_formula = as.formula(paste0(adapted_formula[2], adapted_formula[1], adapted_formula[3]))
    return(adapted_formula)
}

get_model_gene_comparison_analysis_file_path <- function(){
    path = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, 'v_j_comparison_analysis') 
    if (!dir.exists(path)){
        dir.create(path, recursive = true)
    }
    return(path)
}

get_model_gene_comparison_analysis_file_name <- function(){
    path = get_model_gene_comparison_analysis_file_path()
    name = file.path(path, paste0('v_j_comparison_analysis_', MODEL_TYPE, '_', MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND, '_', GENE_WEIGHT_TYPE, '.tsv'))
    return(name)
}

compile_gene_comparison_result <- function(original_model_loglik, gene_var_model_loglik, df_diff, lrt_p, variable_to_vary){
    result = data.table(motif_length_5_end = LEFT_NUC_MOTIF_COUNT, Motif_length_3_end = RIGHT_NUC_MOTIF_COUNT, motif_type = MOTIF_TYPE, gene_weight_type = GENE_WEIGHT_TYPE, upper_bound = UPPER_TRIM_BOUND, lower_bound = LOWER_TRIM_BOUND, model_type = MODEL_TYPE, terminal_melting_5_end_length = LEFT_SIDE_TERMINAL_MELT_LENGTH, original_model_loglik = original_model_loglik, variable_to_vary = variable_to_vary, model_varying_var_by_gene_type_loglik = gene_var_model_loglik, df_diff = df_diff, likelihood_ratio_p = lrt_p) 
    return(result)
}

write_gene_comparison_result_dt <- function(result){
    file_name = get_model_gene_comparison_analysis_file_name() 
    fwrite(result, file_name, sep = '\t')
    return(result)
}
