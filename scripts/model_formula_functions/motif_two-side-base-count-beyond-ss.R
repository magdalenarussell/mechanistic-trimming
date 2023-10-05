stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10

get_start_list <- function(motif_data, trim_type = TRIM_TYPE){
    start_list = NULL 
    return(start_list)
}

source(paste0(MOD_PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_base_count.R'))

get_parameter_vector <- function(trims, genes){
    motif_positions = c()
    for (i in seq(length(trims))){
        motif_positions = c(motif_positions, get_positions(trims[i]))
    }

    left_vars = c()
    right_vars = c()
    for (i in seq(length(trims))){
        left_vars = c(left_vars, get_all_base_variables('5end', trims[i]))
        right_vars = c(right_vars, get_all_base_variables('3end', trims[i]))
    }

    return(c(motif_positions, left_vars, right_vars))
}

get_model_formula <- function(trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    vars = get_parameter_vector(trims, genes)

    base_vars = vars[vars %like% 'base_count']
    base_vars_collapse = paste(base_vars, collapse = ' + ')

    motif_vars = vars[vars %like% 'motif']
    motif_positions_together = paste(motif_vars, collapse = ' + ')

    gene_groups = paste(paste0(genes, '_group'), collapse = ' ,')

    formula = formula(paste0('cbind(weighted_observation, interaction(', gene_groups, ')) ~ ', base_vars_collapse, ' + ', motif_positions_together))

    return(formula)
}

process_single_data_for_model_fit <- function(group_motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    together = process_for_two_side_base_count(group_motif_data, beyond_motif = TRUE, single_stranded = TRUE, gene_type = gene_type, trim_type = trim_type)
    return(together)
}
