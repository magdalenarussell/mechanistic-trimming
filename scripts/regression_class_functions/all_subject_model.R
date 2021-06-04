get_regression_formula <- function(){
    # get grouped column names for all motifs
    motif_positions = c(paste0('motif_pos5_', c(seq(LEFT_NUC_MOTIF_COUNT, 1))), paste0('motif_pos3_', c(seq(1, RIGHT_NUC_MOTIF_COUNT))))
    # sum positions
    motif_positions_sum = paste0(motif_positions, collapse = ' + ')
    motif_positions_sum = paste0(motif_positions_sum, ' + log(p_gene)')
    formula = as.formula(paste('log(p_trim_given_gene)', motif_positions_sum, sep = ' ~ '))
    return(formula)
}

#TODO see notes!!
set_data_factors <- function(data, ref_base){
    motif_positions = c(paste0('pos5_', c(seq(LEFT_NUC_MOTIF_COUNT, 1))), paste0('pos3_', c(seq(1, RIGHT_NUC_MOTIF_COUNT))))
    for (motif_position in motif_positions){
        data[[motif_position]] = as.factor(data[[motif_position]])
        data[[motif_position]] = relevel(data[[motif_position]], ref = ref_base)
    }
    data$gene = as.factor(data$gene)
    return(data)
}


fit_model <- function(data, ref_base = 'A'){
    formula = get_regression_formula()
    # options(contrasts=c('contr.sum','contr.poly'))
    # data = set_data_factors(data, ref_base)
    # model = speedlm(formula, data, sparse = FALSE)
    model = lm(formula, data, weights = p_gene)
    return(model)
}

predict_trimming_dist_given_gene <- function(data, gene_name){
    # filter out unobserved, for now
    data = data[observed == TRUE]

    model = fit_model(data)
    gene_possible_motifs = get_all_motifs_by_gene(gene_name)

    # get population p_gene
    rep_sizes = unique(data[, c('subject', 'count_subject')])
    gene_rep_sizes = unique(data[gene == gene_name][, c('count_subject_gene', 'subject')])
    total_count = sum(rep_sizes$count_subject)
    total_gene_count = sum(gene_rep_sizes$count_subject_gene)
    p_gene = total_gene_count/total_count
    gene_possible_motifs$total_count = total_count 
    gene_possible_motifs$total_gene_count = total_gene_count
    gene_possible_motifs$p_gene = p_gene

    gene_possible_motifs$predicted_prob = exp(predict(model, gene_possible_motifs))

    return(gene_possible_motifs)
}

get_predicted_dist_file_name <- function(){
    path = file.path(OUTPUT_PATH, 'predicted_trimming_distributions')
    dir.create(path)

    complete_path = file.path(path, paste0(GENE_NAME, '_predicted_trimming_dist.tsv'))
    return(complete_path)
}


write_predicted_trimming_dist_by_genes <- function(data, genes){
    all_genes_predicted_data = data.table()
    for (gene in genes){
        temp_file = predict_trimming_dist_given_gene(data, gene)
        all_genes_predicted_data = rbind(all_genes_predicted_data, temp_file)
        print(paste0('finished predictions for ', gene))
    }
    file_name = get_predicted_dist_file_name()
    fwrite(all_genes_predicted_data, file_name, sep = '\t')
}

#get_fitted_PWM_path <- function(){
#    name = paste0(paste0(TRIM_TYPE, '_PWM_',REGRESSION_TYPE, '_', MOTIF_TYPE, '_', PFM_TYPE, '.tsv'))
#    path = file.path('_ignore', 'pfms', 'complete', name)
#    return(path)
#}

#get_contrast_mapping <- function(data, variable){
#    contrasts = contrasts(data[[variable]])
#    entry_map = seq(1, ncol(contrasts))
#    names(entry_map) = rownames(contrasts)[1:ncol(contrasts)]
#    return(entry_map)
#}

#get_fitted_PWM <- function(data, write.table = TRUE){
#    #TODO add test that all contrasts are contrast sums
#    options(contrasts=c('contr.sum','contr.poly'))
#    bases = c('A', 'T', 'C', 'G')
#    positions = c(paste0('pos5_', c(seq(LEFT_NUC_MOTIF_COUNT, 1))), paste0('pos3_', c(seq(1, RIGHT_NUC_MOTIF_COUNT))))

#    PWM = matrix(NA, nrow = 4, ncol = (LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT))
#    colnames(PWM) = positions
#    rownames(PWM) = bases
#    for (base in bases){
#        data = set_data_factors(data, ref_base = base) 
#        model = fit_model(data, ref_base = base)
#        entry_map = get_contrast_mapping(data, 'pos3_1')
#        for (position in positions){
#            for (entry in names(entry_map)){
#                if (is.na(PWM[entry, position])){
#                    PWM[entry, position] = coef(model)[[paste0(position, entry_map[entry])]]
#                } else {
#                    matrix_entry = PWM[entry, position]
#                    coeffiecient_entry = coef(model)[[paste0(position, entry_map[entry])]]
#                    stopifnot(isTRUE(all.equal(matrix_entry, coeffiecient_entry)))
#                }
#            }
#        }
#    }
#    if (isTRUE(write.table)){
#        name = get_fitted_PWM_path()
#        write.table(as.data.frame(PWM), name, sep = '\t')
#    }
#    return(PWM)
#}
        

