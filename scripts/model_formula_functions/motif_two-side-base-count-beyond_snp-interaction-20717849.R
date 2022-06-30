stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)


DATA_GROUP <<- 'ungrouped'

get_start_list <- function(motif_data){
    start_list = NULL 
    return(start_list)
}

source(paste0(PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_base_count.R'))

get_model_formula <- function(){
    left_vars = get_all_base_variables('left')
    right_vars = get_all_base_variables('right')

    left_vars_collapse = paste(left_vars, collapse = ' + ')
    right_vars_collapse = paste(right_vars, collapse = ' + ')

    motif_positions = get_positions() 
    motif_positions_together = paste(motif_positions, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ (', left_vars_collapse, ' + ', right_vars_collapse, ' + ', motif_positions_together, ')*snp'))

    return(formula)
}

process_data_for_model_fit <- function(group_motif_data){
    together = process_for_two_side_base_count(group_motif_data, beyond_motif = TRUE)
    snp = 20717849
    snp_genotypes = compile_all_genotypes_snp_list(snp)  
    colnames(snp_genotypes) = c('localID', 'snp')
    #TODO remove this HIP conversion
    gwas_id = fread(ID_MAPPING_FILE)
    snp_genotypes = merge(snp_genotypes, gwas_id)[, -c('scanID', 'localID')]

    together = merge(together, snp_genotypes, by = 'subject')
    together = calculate_subject_gene_weight(together)
    
    return(together)
}
