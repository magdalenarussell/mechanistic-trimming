source(paste0('scripts/model_group_functions/', MODEL_GROUP, '.R'))

# Aggregate all subject data
aggregate_subject_data_by_trim_gene <- function(subject_data){
    subject_data = calculate_subject_gene_weight(subject_data)
    aggregated_subject = subject_data[, .N, by = .(trim_length, gene, motif, gene_type, gene_weight, subject, observed)]
    setnames(aggregated_subject, 'N', 'count')
    aggregated_subject[observed == FALSE, count := 0]
    return(aggregated_subject)
}

split_motif_column_by_motif_position <- function(aggregated_subject_data){
    split_data = aggregated_subject_data %>% separate('motif', c(paste0('motif_5end_pos', seq(LEFT_NUC_MOTIF_COUNT, 1)), paste0('motif_3end_pos', seq(1, RIGHT_NUC_MOTIF_COUNT))), sep = seq(1, LEFT_NUC_MOTIF_COUNT+RIGHT_NUC_MOTIF_COUNT-1))
    
    together = merge(aggregated_subject_data, split_data)
    return(together)
}

aggregate_all_subject_data <- function(directory = get_subject_motif_output_location()){
    #TODO add ability to look at only one subject or group
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    together = foreach(file = files, .combine=rbind) %dopar% {
        file_data = fread(file)
        print(paste(file))
        aggregated_subject = aggregate_subject_data_by_trim_gene(file_data)
        split_motif_column_by_motif_position(aggregated_subject) 
    }
    stopImplicitCluster()
    return(together)
}
