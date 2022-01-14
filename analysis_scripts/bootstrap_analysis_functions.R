get_bootstrap_sample_from_gene_list <- function(motif_data, genes){
    boot_data = data.table()
    count = nrow(genes)
    for (sample in seq(count)){
        sample_gene = genes[sample]$gene
        temp = motif_data[gene == sample_gene]
        temp$sample_draw = genes[sample]$sample_draw
        boot_data = rbind(boot_data, temp)
    }
    return(boot_data)
}

generate_bootstrap_sample <- function(motif_data, sample_size){
    sample_genes = sample(unique(motif_data$gene), sample_size, replace = TRUE)
    sample_data = data.table()
    sample_draw = 1
    for (sample in sample_genes){
        temp = motif_data[gene == sample]
        temp$sample_draw = sample_draw
        sample_data = rbind(sample_data, temp)
        cols = colnames(sample_data)
        sample_draw = sample_draw + 1
    }
    return(sample_data)
}

edit_formula <- function(formula){
    new_left_side = 'cbind(weighted_observation, interaction(gene, subject, sample_draw)) ~ '
    old_right_side = str_split(formula, '~')[[3]]
    together = paste(new_left_side, old_right_side, collapse = '')
    return(as.formula(together))
}

get_pwm_matrix_file_path_bootstrap <- function(sample_number){
    if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(OUTPUT_PATH, 'bootstrap_genes', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, paste0(MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0('bootstrap_sample_', sample_number), paste0(MODEL_GROUP, '_predicted_coefficient_matrix'), GENE_WEIGHT_TYPE, model)
    dir.create(path, recursive = TRUE)
    
    return(path)
}

get_data_genes_path_bootstrap <- function(sample_number){
    path = file.path(OUTPUT_PATH, 'bootstrap_genes', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, paste0(MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0('bootstrap_sample_', sample_number))
    dir.create(path, recursive = TRUE)
    
    return(path)
}

fit_model_bootstrap <- function(group_motif_data){
    group_motif_data = process_data_for_model_fit(group_motif_data)
    formula = get_model_formula()
    formula =  edit_formula(formula)
    group_motif_data = set_contrasts(group_motif_data)
    stopifnot(MODEL_TYPE == 'motif_distance_two_side_terminal_melting')
    start_list = get_start_list()
    distance_starts = rep(0, UPPER_TRIM_BOUND-LOWER_TRIM_BOUND)
    start = c(start_list, distance_starts, 0, 0)
    model = mclogit(formula, 
                    data = group_motif_data,
                    start = start) 

    return(model)
}

get_coeffiecient_matrix_bootstrap <- function(group_motif_data, ref_base){
    # group_motif_data = process_data_for_model_fit(group_motif_data)
    group_motif_data = set_contrasts(group_motif_data, ref_base)
    model = fit_model_bootstrap(group_motif_data)
    positions = get_positions()

    together = matrix(0, nrow = 4, ncol = LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT)
    colnames(together) = positions
    rownames(together) = c('A', 'C', 'T', 'G')
 
    for (position in positions){
        levels = get_levels(group_motif_data, ref_base, position)
        indices = levels[base %in% c('A', 'C', 'G', 'T') & !is.na(number)]$number
        for (index in indices){
            base = levels[number == index]$base
            together[base, position] = coef(model)[[paste0(position, index)]]
        }
    }

    for (position in positions){
        levels = get_levels(group_motif_data, ref_base, position)
        missing_base = levels[is.na(number)]$base
        together[missing_base,position] = -1*sum(together[, position])
    }

    return(together)
}


fit_model_bootstrap_genes <- function(motif_data, sample_number, write_coeffs = TRUE){
    motif_data = process_data_for_model_fit(motif_data)
    stopifnot('sample_draw' %in% colnames(motif_data))
    gene_path = get_data_genes_path_bootstrap(sample_number)
    genes = motif_data[, .N, by = .(gene, sample_draw)][, -c('N')]
    fwrite(genes, paste0(gene_path, '/sampled_genes.tsv'), sep = '\t')

    if (grepl('motif', MODEL_TYPE, fixed = TRUE)){
        # calculate coeffiecients
        pwm_matrix = get_coeffiecient_matrix_bootstrap(motif_data, ref_base = 'A')
        pwm_dt = as.data.table(pwm_matrix)
        pwm_dt$base = rownames(pwm_matrix)
        pwm_dt$model_group = MODEL_GROUP

        # save coefficients
        pwm_file_path = get_pwm_matrix_file_path_bootstrap(sample_number)
        pwm_file_name = get_pwm_matrix_file_name(subgroup = NULL)
        location = file.path(pwm_file_path, pwm_file_name)

        if (isTRUE(write_coeffs)){
            fwrite(pwm_dt, location, sep = '\t')
        } else {
            return(pwm_dt)
        }
    }
}

