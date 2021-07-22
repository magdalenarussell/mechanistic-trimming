source(paste0('scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'))
source(paste0('scripts/model_group_functions/', MODEL_GROUP, '.R'))

get_positions <- function(){
    positions = c(paste0('motif_5end_pos', seq(LEFT_NUC_MOTIF_COUNT, 1)), paste0('motif_3end_pos', seq(1, RIGHT_NUC_MOTIF_COUNT)))
    return(positions)
}

# Aggregate all subject data
aggregate_subject_data_by_trim_gene <- function(subject_data){
    aggregated_subject = subject_data[, .N, by = .(trim_length, gene, motif, gene_type, subject, observed)]
    setnames(aggregated_subject, 'N', 'count')
    aggregated_subject[observed == FALSE, count := 0]
    return(aggregated_subject)
}

split_motif_column_by_motif_position <- function(aggregated_subject_data){
    positions = get_positions()
    split_data = aggregated_subject_data %>% separate('motif', positions, sep = seq(1, LEFT_NUC_MOTIF_COUNT+RIGHT_NUC_MOTIF_COUNT-1))
    
    together = merge(aggregated_subject_data, split_data)
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

get_model_formula <- function(){
    motif_positions = get_positions() 
    motif_positions_together = paste(motif_positions, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', motif_positions_together))
    return(formula)
}

get_contrasts_list <- function(){
    contrast_list = list()
    positions = get_positions() 

    for (position in positions){
        contrast_list[[position]] = 'contr.sum'
    }
    return(contrast_list)
}

get_start_list <- function(){
    positions = get_positions()
    position_count = length(positions)
    start_list = rep(0, position_count * 3)
    return(start_list)
}

fit_model <- function(group_motif_data){
    formula = get_model_formula()
    contrast_list = get_contrasts_list()
    start_list = get_start_list()

    model = mclogit(formula, 
                    data = group_motif_data, 
                    start = start_list, 
                    contrasts = contrast_list) 

    return(model)
}

get_predicted_dist_file_path <- function(){
    path = file.path(OUTPUT_PATH, 'predicted_trimming_distributions', paste0(MODEL_GROUP,'_', MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT), GENE_WEIGHT_TYPE)
    dir.create(path, recursive = TRUE)
    return(path)
}

relevel_contrasts <- function(group_motif_data, ref_base){
    positions = get_positions()
    for (position in positions){
        group_motif_data[[position]] = as.factor(group_motif_data[[position]])
        group_motif_data[[position]] = relevel(group_motif_data[[position]], ref_base)
    }
    return(group_motif_data)
}

get_levels <- function(group_motif_data, ref_base){
    positions = get_positions()
    options(contrasts = rep("contr.sum", 2))
    group_motif_data = relevel_contrasts(group_motif_data, ref_base)
    contrasts = contrasts(group_motif_data[[positions[1]]])
    levels = data.table(base = rownames(contrasts(group_motif_data[[positions[1]]])), number = c(1, 2, 3, NA))
    return(levels)
}

get_coeffiecient_matrix <- function(group_motif_data, ref_base){
    group_motif_data = relevel_contrasts(group_motif_data, ref_base)
    model = fit_model(group_motif_data)

    levels = get_levels(group_motif_data, ref_base)

    together = matrix(0, nrow = 4, ncol = LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT)
    colnames(together) = get_positions() 
    rownames(together) = c('A', 'C', 'T', 'G')

    for (coef in seq(1, length(coef(model)))){
        name = names(coef(model)[coef])
        position = substring(name, 1, 15)
        num = substring(name, 16, 16)
        base = levels[number == num]$base
        together[base, position] = coef(model)[coef]
    }

    missing_base = levels[is.na(number)]$base
    together[missing_base,] = -1*colSums(together)

    return(together)
}

get_pwm_matrix_file_path <- function(){
    path = file.path(OUTPUT_PATH, 'predicted_coefficient_matrix', paste0(MODEL_GROUP, '_', MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT), GENE_WEIGHT_TYPE)
    dir.create(path, recursive = TRUE)
    
    return(path)
}


