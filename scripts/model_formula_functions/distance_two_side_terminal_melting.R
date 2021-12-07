get_model_formula <- function(){
    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ left_terminal_melting + left_terminal_melting + as.factor(trim_length)'))
    return(formula)
}

get_melting_temp <- function(){
    whole_nucseq = fread(get(paste0('WHOLE_NUCSEQS_', ANNOTATION_TYPE)))
    setnames(whole_nucseq, 'gene', 'gene_names', skip_absent = TRUE)
    setnames(whole_nucseq, 'sequence', 'sequences', skip_absent = TRUE)

    trims = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)
    
    genes = whole_nucseq$gene_names[substring(whole_nucseq$gene_names, 4, 4) == toupper(substring(GENE_NAME, 1, 1))]
    together = data.table(gene = rep(genes, length(trims)), trim_length = rep(trims, length(genes)))
    together = merge(together, whole_nucseq, by.x = 'gene', by.y = 'gene_names')

    setnames(together, 'gene', GENE_NAME)
    map = get_common_genes_from_seqs(together)
    together = merge(together, map, by = GENE_NAME)[,-c('v_gene')]

    if ('name' %in% colnames(together)){
        together = together[, -c('name')]
    }

    #NOT INCLUDING PNUCS IN THIS TERMINAL GC CONTENT CALCULATION
    # get terminal seq
    together[, right_seq := substring(sequences, nchar(sequences) - trim_length + 1, nchar(sequences))]
    together[, left_seq := substring(sequences, nchar(sequences) - (trim_length + LEFT_SIDE_TERMINAL_MELT_LENGTH) + 1, nchar(sequences)-trim_length)]


    # get GC content
    left_seq_list = DNAStringSet(together$left_seq)
    right_seq_list = DNAStringSet(together$right_seq)

    left_base_counts = as.data.table(letterFrequency(left_seq_list, letters="ACGT", OR = 0))
    left_base_counts[(C+G+A+T) < 14 , left_terminal_melting := 4*(C+G)+2*(A+T)]
    left_base_counts[(C+G+A+T) >= 14 , left_terminal_melting := 64.9 + 41*(G+C-16.4)/(A+T+G+C)]
    
    right_base_counts = as.data.table(letterFrequency(right_seq_list, letters="ACGT", OR = 0))
    right_base_counts[(C+G+A+T) < 14 , right_terminal_melting := 4*(C+G)+2*(A+T)]
    right_base_counts[(C+G+A+T) >= 14 , right_terminal_melting := 64.9 + 41*(G+C-16.4)/(A+T+G+C)]
 
    # merge
    together = cbind(together, left_base_counts, right_base_counts)
    return(unique(together[, c('gene', 'trim_length', 'left_terminal_melting', 'right_terminal_melting')]))
}

process_data_for_model_fit <- function(group_motif_data){
    row_count = nrow(group_motif_data)
    if (!('left_terminal_melting' %in% colnames(group_motif_data))){
        terminal_gc = get_melting_temp()
        if (is.factor(group_motif_data$trim_length)){ 
            terminal_gc$trim_length = as.factor(terminal_gc$trim_length)
        }
        together = merge(group_motif_data, terminal_gc, by = c('gene', 'trim_length'))
    } else {
        together = group_motif_data
    }
    stopifnot(nrow(together) == row_count)
    return(together)
}
