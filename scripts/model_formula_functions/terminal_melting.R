get_model_formula <- function(){
    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ terminal_melting'))
    return(formula)
}

get_melting_temp <- function(){
    whole_nucseq = fread(get(paste0('WHOLE_NUCSEQS_', ANNOTATION_TYPE)))
    trims = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)
    
    genes = whole_nucseq$gene_names[substring(whole_nucseq$gene_names, 4, 4) == toupper(substring(GENE_NAME, 1, 1))]
    together = data.table(gene = rep(genes, length(trims)), trim_length = rep(trims, length(genes)))
    together = merge(together, whole_nucseq, by.x = 'gene', by.y = 'gene_names')

    setnames(together, 'gene', GENE_NAME)
    map = get_common_genes_from_seqs(together)
    together = merge(together, map, by = GENE_NAME)[,-c('v_gene')]

    # get length of terminal seq
    #NOT INCLUDING PNUCS IN THIS TERMINAL GC CONTENT CALCULATION
    together[, depth := trim_length + LEFT_NUC_MOTIF_COUNT]

    # get terminal seq
    together[, terminal_seq := substring(sequences, nchar(sequences) - depth + 1, nchar(sequences))]

    # get GC content
    seq_list = DNAStringSet(together$terminal_seq)
    base_counts = as.data.table(letterFrequency(seq_list, letters="ACGT", OR = 0))
    base_counts[(C+G+A+T) < 14 , terminal_melting := 4*(C+G)+2*(A+T)]
    base_counts[(C+G+A+T) >= 14 , terminal_melting := 64.9 + 41*(G+C-16.4)/(A+T+G+C)]
    
    # merge
    together = cbind(together, base_counts)
    return(unique(together[, c('gene', 'trim_length', 'terminal_melting')]))
}

process_data_for_model_fit <- function(group_motif_data){
    row_count = nrow(group_motif_data)
    if (!('terminal_melting' %in% colnames(group_motif_data))){
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

