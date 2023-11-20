stopifnot(LEFT_NUC_MOTIF_COUNT > 0 | RIGHT_NUC_MOTIF_COUNT > 0)
stopifnot(TRIM_TYPE == 'v-j_trim')

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10

get_start_list <- function(motif_data, trim_type = TRIM_TYPE){
    start_list = NULL 
    return(start_list)
}

source(paste0(MOD_PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/two_side_base_count.R'))
source(paste0(MOD_PROJECT_PATH, '/scripts/model_formula_functions/model_formula_specific_functions/mh.R'))

get_parameter_vector <- function(trims, genes){
    motif_positions = c()
    for (i in seq(length(trims))){
        motif_positions = c(motif_positions, get_positions(trims[i]))
    }

    mh_vars = get_all_mh_prop_variables(overlap_vector = c(0, 1, 2, 3, 4))

    left_vars = c()
    right_vars = c()
    for (i in seq(length(trims))){
        left_vars = c(left_vars, paste0(get_all_base_variables('5end', trims[i]), '_prop'))
        right_vars = c(right_vars, paste0(get_all_base_variables('3end', trims[i]), '_prop'))
    }

    return(c(motif_positions, left_vars, right_vars, mh_vars))
}

get_model_formula <- function(trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    vars = get_parameter_vector(trims, genes)

    trim_vars = paste(trims, collapse = ' + ')

    mh_vars = vars[vars %like% 'mh_prop']
    mh_vars_collapse = paste(mh_vars, collapse = ' + ')

    base_vars = vars[vars %like% 'base_count']
    base_vars_collapse = paste(base_vars, collapse = ' + ')

    motif_vars = vars[vars %like% 'motif']
    motif_positions_together = paste(motif_vars, collapse = ' + ')

    gene_groups = paste(paste0(genes, '_group'), collapse = ' ,')

    formula = formula(paste0('cbind(weighted_observation, interaction(', gene_groups, ')) ~ ', base_vars_collapse, ' + ', motif_positions_together, ' + ', mh_vars_collapse, ' + ', trim_vars))

    return(formula)
}

process_single_data_for_model_fit <- function(group_motif_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    together = process_for_two_side_base_count(group_motif_data, beyond_motif = TRUE, whole_nucseq = whole_nucseq, gene_type = gene_type, trim_type = trim_type)
    together = process_for_mh(together, whole_nucseq = whole_nucseq, overlap_vector = c(0, 1, 2, 3, 4), trim_type = trim_type, gene_type = gene_type)
    base_count_cols = colnames(together)[colnames(together) %like% 'base_count']
    base_count_cols = base_count_cols[base_count_cols %like% trim_type]
    base_count_col_types = unique(substring(base_count_cols, 1, 22))

    for (col in base_count_col_types){
        gc_temp = paste0(col, '_GC')
        at_temp = paste0(col, '_AT')
        gc_prop_temp = paste0(col, '_GC_prop')
        at_prop_temp = paste0(col, '_AT_prop')
        together[((get(gc_temp) + get(at_temp)) > 0), paste(gc_prop_temp) := get(gc_temp)/(get(gc_temp) + get(at_temp))]
        together[((get(gc_temp) + get(at_temp)) > 0), paste(at_prop_temp) := get(at_temp)/(get(gc_temp) + get(at_temp))]
        together[((get(gc_temp) + get(at_temp)) == 0), paste(at_prop_temp) := 0]
        cols = colnames(together)[!(colnames(together) %in% c(gc_temp, at_temp))]
        together = together[, ..cols]
    }

    return(together)
}
