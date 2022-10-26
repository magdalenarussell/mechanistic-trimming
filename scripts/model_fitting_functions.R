source(paste0('scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'))
source(paste0('scripts/model_group_functions/', MODEL_GROUP, '.R'))

get_positions <- function(){
    if (LEFT_NUC_MOTIF_COUNT > 0){
        left = paste0('motif_5end_pos', seq(LEFT_NUC_MOTIF_COUNT, 1))    
    } else {
        left = c()
    } 
    if (RIGHT_NUC_MOTIF_COUNT > 0){
        right = paste0('motif_3end_pos', seq(1, RIGHT_NUC_MOTIF_COUNT))
    } else {
        right = c()
    }
    positions = c(left, right)
    return(positions)
}

split_motif_column_by_motif_position <- function(aggregated_subject_data){
    positions = get_positions()

    if (LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT > 0){
        split_data = aggregated_subject_data %>% separate('motif', positions, sep = seq(1, LEFT_NUC_MOTIF_COUNT+RIGHT_NUC_MOTIF_COUNT-1))
        together = merge(aggregated_subject_data, split_data, by = colnames(aggregated_subject_data)[!colnames(aggregated_subject_data) == 'motif'])
    } else {
        together = aggregated_subject_data
    }
    return(together)
}

aggregate_all_subject_data <- function(directory = get_subject_motif_output_location()){
    stopifnot(LEFT_NUC_MOTIF_COUNT <= 10)
    stopifnot(LEFT_SIDE_TERMINAL_MELT_LENGTH <= 10 | is.na(LEFT_SIDE_TERMINAL_MELT_LENGTH))
    desired_file_count = length(list.files(get(paste0('TCR_REPERTOIRE_DATA_', ANNOTATION_TYPE))))
    if (!dir.exists(directory) | !(length(list.files(directory)) == desired_file_count)) {
        print('compiling motif data, first')
        compile_all_data(get(paste0('TCR_REPERTOIRE_DATA_', ANNOTATION_TYPE))) 
    }  
    
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    together = foreach(file = files, .combine=rbind) %dopar% {
        file_data = fread(file)
        print(paste(file))
        file_data
    }
    
    if (MODEL_TYPE %like% 'dna_shape') {
        together = convert_data_to_motifs(together, left_window_size = LEFT_NUC_MOTIF_COUNT + 2, right_window_size = RIGHT_NUC_MOTIF_COUNT + 2)
        processed_motif_data = process_data_for_model_fit(together)
        motif_data = convert_data_to_motifs(processed_motif_data)
    } else {
        together = convert_data_to_motifs(together)
        motif_data = process_data_for_model_fit(together)
    }

    motif_data = motif_data[, -c('left_nucs', 'right_nucs', 'left_motif', 'right_motif')]
    together_pos = split_motif_column_by_motif_position(motif_data) 
    weighted_together = calculate_subject_gene_weight(together_pos)
    stopImplicitCluster()
    return(weighted_together)
}

fit_model <- function(group_motif_data, formula = get_model_formula()){
    stopifnot(unique(group_motif_data$gene_weight_type) == GENE_WEIGHT_TYPE)
    group_motif_data = set_contrasts(group_motif_data)
    start_list = get_start_list(group_motif_data)

    model = mclogit(formula, 
                    data = group_motif_data,
                    start = start_list) 

    return(model)
}

get_predicted_dist_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(MODEL_GROUP, '_predicted_trimming_distributions'), GENE_WEIGHT_TYPE, paste0('motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT), model)
    dir.create(path, recursive = TRUE)
    return(path)
}

get_levels <- function(group_motif_data, ref_base, position){
    group_motif_data = set_contrasts(group_motif_data, ref_base)
    contrasts = contrasts(group_motif_data[[position]])
    levels = data.table(base = rownames(contrasts(group_motif_data[[position]])), number = c(seq(ncol(contrasts(group_motif_data[[position]]))), NA))
    return(levels)
}


get_coeffiecient_matrix <- function(group_motif_data, ref_base, formula = get_model_formula()){
    stopifnot(MODEL_TYPE %like% "motif")
    group_motif_data = set_contrasts(group_motif_data, ref_base)
    model = fit_model(group_motif_data, formula = formula)
    positions = get_positions()

    together = matrix(0, nrow = 4, ncol = LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT)
    if (MODEL_TYPE %like% 'snp-interaction'){
        snp_together = matrix(0, nrow = 4, ncol = LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT)
        colnames(snp_together) = paste0(positions)
        rownames(snp_together) = c('A', 'C', 'T', 'G')
    } else {
        snp_together = NULL
    }

    colnames(together) = positions
    rownames(together) = c('A', 'C', 'T', 'G')
 
    for (position in positions){
        levels = get_levels(group_motif_data, ref_base, position)
        indices = levels[base %in% c('A', 'C', 'G', 'T') & !is.na(number)]$number
        for (index in indices){
            base = levels[number == index]$base
            together[base, position] = coef(model)[[paste0(position, index)]]
            if (MODEL_TYPE %like% 'snp-interaction'){
                snp_together[base, position] = coef(model)[[paste0(position, index, ':snp')]]
            }
        }
    }

    for (position in positions){
        levels = get_levels(group_motif_data, ref_base, position)
        missing_base = levels[is.na(number)]$base
        together[missing_base,position] = -1*sum(together[, position])
        if (MODEL_TYPE %like% 'snp-interaction'){
            snp_together[missing_base,position] = -1*sum(snp_together[, position])
        }
    }

    return(list(result = together, model = model, snp_interaction_result = snp_together))
}

get_complete_distance_coefficients <- function(model_coefficients_dt){
    parameter_names = paste0('trim_length_', seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)) 
    contr_sum_var = parameter_names[length(parameter_names)] 
    parameter_names_subset = parameter_names[!(parameter_names == contr_sum_var)]
    distance_coefs = model_coefficients_dt[parameter %like% 'trim_length']
    distance_coefs$parameter = parameter_names_subset
    missing = data.table(parameter = contr_sum_var, coefficient = -1 * sum(distance_coefs$coefficient), base = NA)
    together = rbind(distance_coefs, missing)
    return(together)
}

format_model_coefficient_output <- function(model, formatted_pwm_matrix = NULL){
    coef_dt = as.data.frame(model$coefficients)
    colnames(coef_dt) = c('coefficient')
    coef_dt$parameter = rownames(coef_dt)
    coef_dt = as.data.table(coef_dt)
    if (!is.null(formatted_pwm_matrix)){
        stopifnot(MODEL_TYPE %like% 'motif')
        coef_dt = coef_dt[!(parameter %like% 'motif')]
        positions = get_positions()
        not_cols = colnames(formatted_pwm_matrix)[!(colnames(formatted_pwm_matrix)%in% positions)]
        # formatted_pwm = formatted_pwm_matrix[, -c('model_group')] %>% 
        formatted_pwm = formatted_pwm_matrix %>% 
            pivot_longer(!not_cols, values_to = 'coefficient', names_to = 'parameter') %>%
            as.data.table()
        coef_dt = rbind(coef_dt, formatted_pwm, fill = TRUE)
    } else {
        if (MODEL_TYPE %like% 'base-count'){
            coef_dt[parameter %like% 'base_count', base := substring(parameter, nchar(parameter)-1, nchar(parameter))]
            coef_dt[parameter %like% 'base_count', parameter := str_remove(parameter, substring(parameter, nchar(parameter)-2, nchar(parameter)))]
        } else {
            coef_dt$base = NA
        }
    }

    if ((MODEL_TYPE %like% 'distance') & !(MODEL_TYPE %like% 'linear-distance')) {
        dist_coefs = get_complete_distance_coefficients(coef_dt)
        coef_dt = coef_dt[!(parameter %like% 'trim_length')]
        coef_dt = rbind(coef_dt, dist_coefs)
    }
    return(coef_dt)
}

get_pwm_matrix_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(MODEL_GROUP, '_predicted_coefficient_matrix'), GENE_WEIGHT_TYPE, paste0('motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT), model)
    dir.create(path, recursive = TRUE)
    
    return(path)
}

get_coefficient_output_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(MODEL_GROUP, '_predicted_coefficients'), GENE_WEIGHT_TYPE, paste0('motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT), model)
    dir.create(path, recursive = TRUE)
    
    return(path)
}

cluster_sample <- function(motif_data){
    # sample size is the number of gene, subject combos
    motif_data$cluster = interaction(motif_data$subject, motif_data$gene)
    sample_genes = sample(unique(motif_data$cluster), length(unique(motif_data$cluster)), replace = TRUE)
    sample_genes_dt = data.table(sample_genes)
    counts = sample_genes_dt[, .N, by = sample_genes]
    motif_data = merge(motif_data, counts, by.x = 'cluster', by.y = 'sample_genes')
    # updating the total_tcr, p_gene, gene_weight_type, and weighted_observation variables for the newly sampled datasets
    source(paste0('scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'), local = TRUE)
    sample_data = calculate_subject_gene_weight(motif_data)
    return(sample_data)
}

cluster_bootstrap_model_fit <- function(motif_data, formula = get_model_formula(), iter){
    set.seed(66)
    results = data.table()
    for (i in seq(iter)){
        sample_data = cluster_sample(motif_data)
        if (MODEL_TYPE %like% 'motif'){
            result = get_coeffiecient_matrix(sample_data, formula = formula, ref_base = 'A')
            pwm_matrix = result$result
            pwm_dt = as.data.table(pwm_matrix)
            pwm_dt$base = rownames(pwm_matrix)
            if (MODEL_TYPE %like% 'snp-interaction'){
                pwm_snp_matrix = result$snp_interaction_result
                pwm_snp_dt = as.data.table(pwm_snp_matrix)
                pwm_snp_dt$base = rownames(pwm_snp_matrix)
                pwm_snp_dt$snp_interaction = TRUE
                pwm_dt$snp_interaction = FALSE
                pwm_dt = rbind(pwm_dt, pwm_snp_dt)
            } else {
                pwm_dt$snp_interaction = FALSE
            }
            model = result$model
        } else {
            pwm_dt = NULL
            model = fit_model(sample_data, formula)
        }
        coefs = format_model_coefficient_output(model, pwm_dt)
        if (!('snp_interaction' %in% colnames(coefs))){
            coefs$snp_interaction = FALSE
        }
        coefs$iteration = i
        results = rbind(results, coefs)
    }
    results[snp_interaction == TRUE, parameter := paste0(parameter, ':snp')]
    return(results)
}

get_coef_pvalues <- function(bootstrap_results, original_model_results){
    sd_coeff = bootstrap_results[, sd(coefficient), by = .(parameter, base)]
    setnames(sd_coeff, 'V1', 'sd')

    together = merge(original_model_results, sd_coeff)

    together[, zstat := coefficient/sd]
    together[, pvalue := 2*pnorm(-abs(zstat))]
    together$iterations = max(bootstrap_results$iteration)
    return(together)
}
get_model_bootstrap_path <- function(){
    path = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, 'model_bootstrap') 
    if (!dir.exists(path)){
        dir.create(path, recursive = TRUE)
    }
    return(path)
}

get_model_bootstrap_file_name <- function(){
    path = get_model_bootstrap_path()
    name = file.path(path, paste0('bootstrap_', MODEL_TYPE, '_', MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND, '_', GENE_WEIGHT_TYPE, '.tsv'))
    return(name)
}

get_model_name <- function(){
    dir = file.path(PROJECT_PATH, 'models', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY) 
    dir.create(dir)
    file = file.path(dir, paste0(MODEL_TYPE, '_', MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND, '_', GENE_WEIGHT_TYPE, '.rds'))
    return(file)
}

save_model <- function(model){
    file = get_model_name()
    saveRDS(model, file)
}

load_model <- function(){
    file = get_model_name()
    model = readRDS(file)
    return(model)
}
