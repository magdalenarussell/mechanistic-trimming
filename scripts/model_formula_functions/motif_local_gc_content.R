get_model_formula <- function(){
    motif_positions = get_positions() 
    motif_positions_together = paste(motif_positions, collapse = ' + ')

    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ ', motif_positions_together, ' + motif_gc_content'))
    return(formula)
}

get_GC_content <- function(seq_list){
    seq_list = DNAStringSet(seq_list)
    base_counts = as.data.table(letterFrequency(seq_list, letters="ACGT", OR = 0))
    base_counts[, GC_freq := (C+G)/(A+C+G+T)]
    return(base_counts$GC_freq)
}

process_data_for_model_fit <- function(group_motif_data){
    group_motif_data$motif_gc_content = get_GC_content(group_motif_data$motif)
    return(group_motif_data)
}
