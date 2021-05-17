get_regression_formula <- function(conditioning = 'gene'){
    all_motif_sums = c()
    for (motif_col in c('motif', paste0('motif_trim_', seq(2, 18)))){
        # get grouped column names for all motifs
        assign(paste0(motif_col, '_positions'), c(paste0(motif_col, '_pos5_', c(seq(LEFT_NUC_MOTIF_COUNT, 1))), paste0(motif_col, '_pos3_', c(seq(1, RIGHT_NUC_MOTIF_COUNT)))))
        # sum positions
        assign(paste0(motif_col, '_positions_sum'), paste0(get(paste0(motif_col, '_positions')), collapse = ' + '))
        # create weight formula
        assign(paste0(motif_col, '_positions_formula'), paste0('exp(', get(paste0(motif_col, '_positions_sum')), ')'))
        # create vector of all motif formulas
        all_motif_sums = c(all_motif_sums, get(paste0(motif_col, '_positions_formula')))
    }
     
    model_denominator = paste0(all_motif_sums[-1], collapse = ' + ')
    model = paste0('(', motif_positions_sum, ')-log(', model_denominator, ')')
    formula = as.formula(paste('log(p_trim_given_gene)', model, sep = ' ~ '))

    return(formula)
}

#TODO see notes!!
#
set_data_factors <- function(data, ref_base){
    motif_positions = c(paste0('pos5_', c(seq(LEFT_NUC_MOTIF_COUNT, 1))), paste0('pos3_', c(seq(1, RIGHT_NUC_MOTIF_COUNT))))
    for (motif_position in motif_positions){
        data[[motif_position]] = as.factor(data[[motif_position]])
        data[[motif_position]] = relevel(data[[motif_position]], ref = ref_base)
    }
    data$gene = as.factor(data$gene)
    return(data)
}


fit_model <- function(data, conditioning = 'gene', ref_base = 'A'){
    formula = get_regression_formula(conditioning)
    options(contrasts=c('contr.sum','contr.poly'))
    data = set_data_factors(data, ref_base)
    model = speedlm(formula, data, sparse = FALSE)
    return(model)
}

get_fitted_PWM_path <- function(conditioning){
    name = paste0(paste0('PWM_',REGRESSION_TYPE, '_', MOTIF_TYPE, '_', PFM_TYPE, '_', conditioning, '_conditioning.tsv'))
    path = file.path('_ignore', 'pfms', 'complete', name)
    return(path)
}

get_contrast_mapping <- function(data, variable){
    contrasts = contrasts(data[[variable]])
    entry_map = seq(1, ncol(contrasts))
    names(entry_map) = rownames(contrasts)[1:ncol(contrasts)]
    return(entry_map)
}

get_fitted_PWM <- function(data, conditioning = 'gene', write.table = TRUE){
    #TODO add test that all contrasts are contrast sums
    options(contrasts=c('contr.sum','contr.poly'))
    bases = c('A', 'T', 'C', 'G')
    positions = c(paste0('pos5_', c(seq(LEFT_NUC_MOTIF_COUNT, 1))), paste0('pos3_', c(seq(1, RIGHT_NUC_MOTIF_COUNT))))

    PWM = matrix(NA, nrow = 4, ncol = (LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT))
    colnames(PWM) = positions
    rownames(PWM) = bases
    for (base in bases){
        data = set_data_factors(data, ref_base = base) 
        model = fit_model(data, conditioning, ref_base = base)
        entry_map = get_contrast_mapping(data, 'pos3_1')
        for (position in positions){
            for (entry in names(entry_map)){
                if (is.na(PWM[entry, position])){
                    PWM[entry, position] = coef(model)[[paste0(position, entry_map[entry])]]
                } else {
                    matrix_entry = PWM[entry, position]
                    coeffiecient_entry = coef(model)[[paste0(position, entry_map[entry])]]
                    stopifnot(isTRUE(all.equal(matrix_entry, coeffiecient_entry)))
                }
            }
        }
    }
    if (isTRUE(write.table)){
        name = get_fitted_PWM_path(conditioning)
        write.table(as.data.frame(PWM), name, sep = '\t')
    }
    return(PWM)
}
        

