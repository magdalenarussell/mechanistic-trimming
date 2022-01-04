source(paste0('scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'))
source(paste0('scripts/model_group_functions/', MODEL_GROUP, '.R'))
source(paste0('scripts/model_formula_functions/', MODEL_TYPE, '.R'))

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

# Aggregate all subject data
aggregate_subject_data_by_trim_gene <- function(subject_data){
    aggregated_subject = subject_data[, sum(count), by = .(trim_length, gene, motif, gene_type, subject, observed)]
    setnames(aggregated_subject, 'V1', 'count')
    return(aggregated_subject)
}

split_motif_column_by_motif_position <- function(aggregated_subject_data){
    positions = get_positions()

    if (LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT > 0){
        split_data = aggregated_subject_data %>% separate('motif', positions, sep = seq(1, LEFT_NUC_MOTIF_COUNT+RIGHT_NUC_MOTIF_COUNT-1))
        together = merge(aggregated_subject_data, split_data)
    } else {
        together = aggregated_subject_data
    }
    return(together)
}

aggregate_all_subject_data <- function(directory = get_subject_motif_output_location()){
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    together = foreach(file = files, .combine=rbind) %dopar% {
        file_data = fread(file)
        print(paste(file))
        aggregated_subject = aggregate_subject_data_by_trim_gene(file_data)
        split_motif_column_by_motif_position(aggregated_subject) 
    }
    together = calculate_subject_gene_weight(together)
    stopImplicitCluster()
    return(together)
}

get_start_list <- function(){
    positions = get_positions()
    position_count = length(positions)
    start_list = rep(0, position_count * 3)
    return(start_list)
}

fit_model <- function(group_motif_data){
    group_motif_data = process_data_for_model_fit(group_motif_data)
    formula = get_model_formula()
    group_motif_data = set_contrasts(group_motif_data)
    # start_list = get_start_list()

    model = mclogit(formula, 
                    data = group_motif_data) 
                    # start = start_list) 

    return(model)
}

get_predicted_dist_file_path <- function(){
    if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, paste0(MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(MODEL_GROUP, '_predicted_trimming_distributions'), GENE_WEIGHT_TYPE, model)

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
    # group_motif_data = process_data_for_model_fit(group_motif_data)
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

get_pwm_matrix_file_path <- function(){
    if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
        model = paste0(MODEL_TYPE, '_', LEFT_SIDE_TERMINAL_MELT_LENGTH, '_length_melting_left')
    } else {
        model = MODEL_TYPE
    }
    path = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, paste0(MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(MODEL_GROUP, '_predicted_coefficient_matrix'), GENE_WEIGHT_TYPE, model)
    dir.create(path, recursive = TRUE)
    
    return(path)
}
