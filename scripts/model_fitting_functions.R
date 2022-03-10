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

fit_model <- function(group_motif_data){
    stopifnot(unique(group_motif_data$gene_weight_type) == GENE_WEIGHT_TYPE)
    formula = get_model_formula()
    group_motif_data = set_contrasts(group_motif_data)
    start_list = get_start_list(group_motif_data)

    model = mclogit(formula, 
                    data = group_motif_data,
                    start = start_list) 

    return(model)
}

get_predicted_dist_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE)){
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


get_coeffiecient_matrix <- function(group_motif_data, ref_base){
    group_motif_data = set_contrasts(group_motif_data, ref_base)
    model = fit_model(group_motif_data)
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
        formatted_pwm = formatted_pwm_matrix[, -c('model_group')] %>% 
            pivot_longer(!base, values_to = 'coefficient', names_to = 'parameter') %>%
            as.data.table()
        coef_dt = rbind(coef_dt, formatted_pwm, fill = TRUE)
    } else {
        coef_dt$base = NA
    }
    if (MODEL_TYPE %like% 'distance') {
        dist_coefs = get_complete_distance_coefficients(coef_dt)
        coef_dt = coef_dt[!(parameter %like% 'trim_length')]
        coef_dt = rbind(coef_dt, dist_coefs)
    }
    return(coef_dt)
}

get_pwm_matrix_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(MODEL_GROUP, '_predicted_coefficient_matrix'), GENE_WEIGHT_TYPE, paste0('motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT), model)
    dir.create(path, recursive = TRUE)
    
    return(path)
}

get_coefficient_output_file_path <- function(){
    if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(MODEL_GROUP, '_predicted_coefficients'), GENE_WEIGHT_TYPE, paste0('motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT), model)
    dir.create(path, recursive = TRUE)
    
    return(path)
}
