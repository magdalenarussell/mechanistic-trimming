UPPER_TRIM_BOUND <<- 18
LOWER_TRIM_BOUND <<- NA

get_unobserved_motifs <- function(tcr_dataframe){
    tcr_dataframe_observed = tcr_dataframe[,.N, by = .(gene, trim_length)]
    unobserved_df = data.frame()
    for (gene_name in unique(tcr_dataframe_observed$gene)){
        observed = unique(tcr_dataframe_observed[gene == paste(gene_name)]$trim_length)
        unobserved = setdiff(seq(1,20), observed)
        unobserved_df = rbind(unobserved_df, data.frame(gene = rep(gene_name, length(unobserved)), trim_length = unobserved))
    }
    cols = c('whole_seq', 'gene')
    tcr_dataframe = unique(tcr_dataframe[,..cols])
    together = as.data.table(merge(unobserved_df, tcr_dataframe, by = 'gene'))
    together = together[, motif:= get_motif_context_unobserved(whole_seq, trim_length)] 
    together$observed = FALSE
    return(together[,-c('whole_seq')])
}

get_motifs <- function(tcr_dataframe, subject_id){
    tcr_dataframe = tcr_dataframe[get(TRIM_TYPE) > 0]
    cdr3_variable = paste0('cdr3_nucseq_from_', substring(TRIM_TYPE, 1, 1))
    cols = c(paste(cdr3_variable), 'sequences', paste(GENE_NAME), paste(TRIM_TYPE))
    tcr_dataframe = tcr_dataframe[,..cols]
    colnames(tcr_dataframe) = c('trimmed_seq', 'whole_seq', 'gene', 'trim_length')
    tcr_dataframe = tcr_dataframe[,motif:=get_motif_context(whole_seq, trimmed_seq, trim_length)]

    tcr_dataframe$observed = TRUE
    unobserved = get_unobserved_motifs(tcr_dataframe)
    together = rbind(tcr_dataframe[,-c(1,2)], unobserved)
    together$gene_type = GENE_NAME
    together$subject = subject_id
    return(together)
}

set_contrasts <- function(group_motif_data, ref_base = 'A'){
    positions = get_positions()
    for (position in positions){
        group_motif_data[[position]] = as.factor(group_motif_data[[position]])
        group_motif_data[[position]] = relevel(group_motif_data[[position]], ref_base)
        unique_bases = length(unique(group_motif_data[[position]]))
        contrasts(group_motif_data[[position]]) = contr.sum(unique_bases)
        if (unique_bases > 4){
            contrasts = contrasts(group_motif_data[[position]])
            non_base_matrix_pos = which(!(rownames(contrasts) %in% c('A', 'T', 'C', 'G')))
            stopifnot(non_base_matrix_pos < length(rownames(contrasts)))
            contrasts[non_base_matrix_pos,] = rep(0, length(contrasts[non_base_matrix_pos,]))
            contrasts[,non_base_matrix_pos] = rep(0, length(contrasts[,non_base_matrix_pos]))
            contrasts(group_motif_data[[position]]) = contrasts
        }
    }
    return(group_motif_data)
}



