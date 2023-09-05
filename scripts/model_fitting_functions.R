source(paste0(MOD_PROJECT_PATH,'/scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_group_functions/', MODEL_GROUP, '.R'))

get_positions <- function(trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)

    if (LEFT_NUC_MOTIF_COUNT > 0){
        left = c()
        for (i in seq(length(trims))){
            left = c(left, paste0(trims[i], '_motif_5end_pos', seq(LEFT_NUC_MOTIF_COUNT, 1)))
        }
    } else {
        left = c()
    } 
    if (RIGHT_NUC_MOTIF_COUNT > 0){
        right = c()
        for (i in seq(length(trims))){
            right = c(right, paste0(trims[i], '_motif_3end_pos', seq(1, RIGHT_NUC_MOTIF_COUNT)))
        }
    } else {
        right = c()
    }
    positions = c(left, right)
    return(positions)
}

aggregate_all_subject_data <- function(directory = get_subject_motif_output_location(), gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    stopifnot(LEFT_NUC_MOTIF_COUNT <= 10)
    stopifnot(LEFT_SIDE_TERMINAL_MELT_LENGTH <= 10 | is.na(LEFT_SIDE_TERMINAL_MELT_LENGTH))
    desired_file_count = length(list.files(get(paste0('TCR_REPERTOIRE_DATA_', ANNOTATION_TYPE))))
    if (!dir.exists(directory) | !(length(list.files(directory)) == desired_file_count)) {
        print('compiling motif data, first')
        compile_all_data(get(paste0('TCR_REPERTOIRE_DATA_', ANNOTATION_TYPE)), gene_type = gene_type, trim_type = trim_type) 
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
        processed_motif_data = process_data_for_model_fit(together, gene_type = gene_type, trim_type = trim_type)
        motif_data = convert_data_to_motifs(processed_motif_data)
    } else {
        together = convert_data_to_motifs(together)
        motif_data = process_data_for_model_fit(together, gene_type = gene_type, trim_type = trim_type)
    }

    cols = colnames(motif_data)[!(colnames(motif_data) %like% 'left_nucs')]
    cols = cols[!(cols %like% 'right_nucs')]

    motif_data = motif_data[, ..cols]
    together_pos = split_motif_column_by_motif_position(motif_data, trim_type) 
    weighted_together = calculate_subject_gene_weight(together_pos, gene_type = gene_type, trim_type = trim_type)
    stopImplicitCluster()
    return(weighted_together)
}

fit_model <- function(group_motif_data, formula = get_model_formula(TRIM_TYPE, GENE_NAME), trim_type = TRIM_TYPE){
    stopifnot(unique(group_motif_data$gene_weight_type) == GENE_WEIGHT_TYPE)
    group_motif_data = set_contrasts(group_motif_data, trim_type = trim_type)
    start_list = get_start_list(group_motif_data, trim_type)

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
    path = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(MODEL_GROUP, '_predicted_trimming_distributions'), GENE_WEIGHT_TYPE, paste0('motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT), model)
    dir.create(path, recursive = TRUE)
    return(path)
}

get_levels <- function(group_motif_data, ref_base, position, trim_type = TRIM_TYPE){
    group_motif_data = set_contrasts(group_motif_data, ref_base, trim_type)
    contrasts = contrasts(group_motif_data[[position]])
    levels = data.table(base = rownames(contrasts(group_motif_data[[position]])), number = c(seq(ncol(contrasts(group_motif_data[[position]]))), NA))
    return(levels)
}


get_coefficient_matrix <- function(group_motif_data, ref_base, formula = get_model_formula(TRIM_TYPE, GENE_NAME), trim_type = TRIM_TYPE){
    stopifnot(MODEL_TYPE %like% "motif")
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    group_motif_data = set_contrasts(group_motif_data, ref_base, trim_type)
    model = fit_model(group_motif_data, formula = formula, trim_type)
    positions = get_positions(trim_type)

    result = list()
    snp_interaction_result = list()
    for (t in trims){
        subset_pos = positions[positions %like% t]

        together = matrix(0, nrow = 4, ncol = LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT)
        if (MODEL_TYPE %like% 'snp-interaction'){
            snp_together = matrix(0, nrow = 4, ncol = LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT)
            colnames(snp_together) = paste0(subset_pos)
            rownames(snp_together) = c('A', 'C', 'T', 'G')
        } else {
            snp_together = NULL
        }

        colnames(together) = subset_pos 
        rownames(together) = c('A', 'C', 'T', 'G')
 
        for (position in subset_pos){
            levels = get_levels(group_motif_data, ref_base, position, t)
            indices = levels[base %in% c('A', 'C', 'G', 'T') & !is.na(number)]$number
            for (index in indices){
                base = levels[number == index]$base
                together[base, position] = coef(model)[[paste0(position, index)]]
                if (MODEL_TYPE %like% 'snp-interaction'){
                    snp_together[base, position] = coef(model)[[paste0(position, index, ':snp')]]
                }
            }
        }

        for (position in subset_pos){
            levels = get_levels(group_motif_data, ref_base, position, t)
            missing_base = levels[is.na(number)]$base
            together[missing_base,position] = -1*sum(together[, position])
            if (MODEL_TYPE %like% 'snp-interaction'){
                snp_together[missing_base,position] = -1*sum(snp_together[, position])
            }
        }
        result[[t]] = together
        snp_interaction_result[[t]] = snp_together
    }

    return(list(result = result, model = model, snp_interaction_result = snp_interaction_result))
}

#TODO this is not generalized for double trimming
get_complete_distance_coefficients <- function(model_coefficients_dt, trim_type = TRIM_TYPE){
    parameter_names = paste0(trim_type, '_', seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)) 
    contr_sum_var = parameter_names[length(parameter_names)] 
    parameter_names_subset = parameter_names[!(parameter_names == contr_sum_var)]
    distance_coefs = model_coefficients_dt[parameter %like% trim_type]
    distance_coefs$parameter = parameter_names_subset
    missing = data.table(parameter = contr_sum_var, coefficient = -1 * sum(distance_coefs$coefficient), base = NA)
    together = rbind(distance_coefs, missing)
    return(together)
}

format_model_coefficient_output <- function(model, formatted_pwm_matrix = NULL, trim_type = TRIM_TYPE){
    coef_dt = as.data.frame(model$coefficients)
    colnames(coef_dt) = c('coefficient')
    coef_dt$parameter = rownames(coef_dt)
    coef_dt = as.data.table(coef_dt)
    if (!is.null(formatted_pwm_matrix)){
        stopifnot(MODEL_TYPE %like% 'motif')
        coef_dt = coef_dt[!(parameter %like% 'motif')]
        positions = get_positions(trim_type)
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
        dist_coefs = get_complete_distance_coefficients(coef_dt, trim_type = trim_type)
        coef_dt = coef_dt[!(parameter %like% trim_type)]
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
    path = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(MODEL_GROUP, '_predicted_coefficient_matrix'), GENE_WEIGHT_TYPE, paste0('motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT), model)
    dir.create(path, recursive = TRUE)
    
    return(path)
}

get_coefficient_output_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(MODEL_GROUP, '_predicted_coefficients'), GENE_WEIGHT_TYPE, paste0('motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT), model)
    dir.create(path, recursive = TRUE)
    
    return(path)
}

cluster_sample <- function(motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    genes = get_gene_order(gene_type)
    # sample size is the number of gene, subject combos
    if (length(genes) == 1){
        motif_data$cluster = interaction(motif_data$subject, motif_data[[paste0(gene_type, '_group')]])
    } else if (length(genes) == 2){
        motif_data$cluster = interaction(motif_data$subject, motif_data[[paste0(genes[1], '_group')]], motif_data[[paste0(genes[2], '_group')]])
    }

    sample_genes = sample(unique(motif_data$cluster), length(unique(motif_data$cluster)), replace = TRUE)
    sample_genes_dt = data.table(sample_genes)
    counts = sample_genes_dt[, .N, by = sample_genes]
    motif_data = merge(motif_data, counts, by.x = 'cluster', by.y = 'sample_genes')
    # updating the total_tcr, p_gene, gene_weight_type, and weighted_observation variables for the newly sampled datasets
    source(paste0(MOD_PROJECT_PATH,'/scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'), local = TRUE)
    sample_data = calculate_subject_gene_weight(motif_data, gene_type = gene_type, trim_type = trim_type)
    return(sample_data)
}

cluster_bootstrap_model_fit <- function(motif_data, formula = get_model_formula(TRIM_TYPE, GENE_NAME), iter, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    set.seed(66)
    results = data.table()
    for (i in seq(iter)){
        sample_data = cluster_sample(motif_data, gene_type = gene_type, trim_type = trim_type)
        if (MODEL_TYPE %like% 'motif'){
            result = get_coefficient_matrix(sample_data, formula = formula, ref_base = 'A', trim_type = trim_type)
            pwm_matrix = result$result
            pwm_dt = as.data.table(pwm_matrix)
            setnames(pwm_dt, colnames(pwm_dt), map_chr(str_split(colnames(pwm_dt), '\\.', 2), 2))
            if (class(pwm_matrix) == 'list'){
                pwm_dt$base = rownames(pwm_matrix[[1]])
            } else {
                pwm_dt$base = rownames(pwm_matrix)
            }
            if (MODEL_TYPE %like% 'snp-interaction'){
                pwm_snp_matrix = result$snp_interaction_result
                pwm_snp_dt = as.data.table(pwm_snp_matrix)
                if (class(pwm_snp_matrix) == 'list'){
                    pwm_snp_dt$base = rownames(pwm_snp_matrix[[1]])
                } else {
                    pwm_snp_dt$base = rownames(pwm_snp_matrix)
                }
                pwm_snp_dt$snp_interaction = TRUE
                pwm_dt$snp_interaction = FALSE
                pwm_dt = rbind(pwm_dt, pwm_snp_dt)
            } else {
                pwm_dt$snp_interaction = FALSE
            }
            model = result$model
        } else {
            pwm_dt = NULL
            model = fit_model(sample_data, formula, trim_type)
        }
        coefs = format_model_coefficient_output(model, pwm_dt, trim_type)
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

subsample <- function(motif_data, prop, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    stopifnot(TRIM_TYPE %in% c('v_trim', 'j_trim'))
    stopifnot(MODEL_TYPE == 'motif_two-side-base-count-beyond')
    genes = get_gene_order(gene_type)
    trims = get_trim_order(trim_type)

    # sample size is the number of gene, subject combos
    if (length(genes) == 1){
        motif_data$cluster = interaction(motif_data$subject, motif_data[[paste0(gene_type, '_group')]])
    } else if (length(genes) == 2){
        motif_data$cluster = interaction(motif_data$subject, motif_data[[paste0(genes[1], '_group')]], motif_data[[paste0(genes[2], '_group')]])
    }

    # sample proportion of individuals
    size = ceiling(prop * length(unique(motif_data$subject)))
    sample_subj = sample(unique(motif_data$subject), size, replace = FALSE)
    subj_subset_orig = motif_data[subject %in% sample_subj]

    # sample proportion of sequences for each individual
    subj_subset_orig[, subsample_total_tcr := ceiling(total_tcr*prop)]
    cols = c(paste0(genes, '_group'), trims, paste0(trim_type, '_observed'), paste0(trims, '_left_base_count_AT'), paste0(trims, '_left_base_count_GC'), paste0(trims, '_right_base_count_AT'), paste0(trims, '_right_base_count_GC'), paste0(trims, '_motif'), paste0(trims, '_motif_5end_pos1'), paste0(trims, '_motif_3end_pos1'), paste0(trims, '_motif_3end_pos2'), 'subject', 'count', 'subsample_total_tcr', 'cluster')
    subj_subset_orig[, row := seq(1, .N), by = .(subject)]
    subj_subset = subj_subset_orig[subj_subset_orig[, sample(.I, subsample_total_tcr, replace = TRUE, prob = count), by = subject]$V1]
    subj_subset[, count := .N, by = .(subject, row)]
    subj_subset_final = unique(subj_subset[, ..cols])

    # fill in unobserved seq cases
    subj_subset_orig_small = subj_subset_orig[, ..cols][, -c('count')]
    subj_subset_final_small = subj_subset_final[, ..cols][, -c('count')]
    unsampled = fsetdiff(subj_subset_orig_small, subj_subset_final_small) 
    unsampled$count = 0
    subj_subset_final = rbind(subj_subset_final, unsampled)

    # updating the total_tcr, p_gene, gene_weight_type, and weighted_observation variables for the newly sampled datasets
    source(paste0(MOD_PROJECT_PATH,'/scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'), local = TRUE)
    sample_data = calculate_subject_gene_weight(subj_subset_final, gene_type = gene_type, trim_type = trim_type)
    return(sample_data)
}

subsample_model_fit <- function(motif_data, formula = get_model_formula(TRIM_TYPE, GENE_NAME), iter, prop, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    set.seed(66)
    results = data.table()
    for (i in seq(iter)){
        sample_data = subsample(motif_data, prop, trim_type = trim_type, gene_type = GENE_NAME)
        result = get_coefficient_matrix(sample_data, formula = formula, ref_base = 'A', trim_type = trim_type)
        pwm_matrix = result$result
        pwm_dt = as.data.table(pwm_matrix)
        setnames(pwm_dt, colnames(pwm_dt), map_chr(str_split(colnames(pwm_dt), '\\.', 2), 2))
        if (class(pwm_matrix) == 'list'){
            pwm_dt$base = rownames(pwm_matrix[[1]])
        } else {
            pwm_dt$base = rownames(pwm_matrix)
        }
        model = result$model
        coefs = format_model_coefficient_output(model, pwm_dt, trim_type)
        coefs$iteration = i
        results = rbind(results, coefs)
    }
    results$prop = prop
    return(results)
}

get_model_bootstrap_path <- function(){
    path = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE,PARAM_GROUP, 'model_bootstrap') 
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

get_model_subsample_path <- function(){
    path = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, 'model_subsample') 
    if (!dir.exists(path)){
        dir.create(path, recursive = TRUE)
    }
    return(path)
}

get_model_subsample_file_name <- function(prop){
    path = get_model_subsample_path()
    name = file.path(path, paste0('subsample_', prop, '_prop_', MODEL_TYPE, '_', MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND, '_', GENE_WEIGHT_TYPE, '.tsv'))
    return(name)
}

get_model_name <- function(){
    dir = file.path(MOD_PROJECT_PATH, 'models', ANNOTATION_TYPE, PARAM_GROUP) 
    dir.create(dir, recursive = TRUE)
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
