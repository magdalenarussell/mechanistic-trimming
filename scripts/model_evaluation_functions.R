generate_hold_out_sample <- function(motif_data, subject_sample_size){
    set.seed(66)
    sample = sample(unique(motif_data$subject), subject_sample_size, replace = FALSE) 
    sample_data = motif_data[subject == sample] 
    motif_data_subset = motif_data[subject != sample]
    return(list(sample = sample_data, motif_data_subset = motif_data_subset))
}

calculate_cond_log_loss <- function(model, sample_data){
    sample_data$prediction = predict(model, newdata = sample_data, type = 'response')
    sample_data$log_prediction = log(sample_data$prediction)
    sample_data$weighted_log_prediction = sample_data$log_prediction * sample_data$count
    log_loss = -sum(sample_data$weighted_log_prediction)
    return(log_loss)
}

get_model_evaluation_file_name <- function(){
    name = file.path(OUTPUT_PATH, 'model_evaluation.tsv')
    return(name)
}

write_result_dt <- function(log_loss, sample_data){
    file_name = get_model_evaluation_file_name()
    held_out_subject = paste(unique(sample_data$subject), collapse = ', ')
    result = data.table(held_out_subject = held_out_subject, motif_length_5_end = LEFT_NUC_MOTIF_COUNT, motif_length_3_end = RIGHT_NUC_MOTIF_COUNT, log_loss = log_loss, motif_type = MOTIF_TYPE, gene_weight_type = GENE_WEIGHT_TYPE) 
    
    if (file.exists(file_name)){
        results = fread(file_name)
        together = rbind(results, result)
        fwrite(together, file_name, sep = '\t')
    } else {
        fwrite(result, file_name, sep = '\t')
    }
    return(result)
}
