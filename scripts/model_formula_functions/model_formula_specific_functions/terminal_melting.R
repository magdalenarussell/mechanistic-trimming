simple_terminal_melting_calculation <- function(sequence_list){
    # simple terminal melting calculation 
    base_counts = as.data.table(letterFrequency(seq_list, letters="ACGT", OR = 0))
    base_counts[(C+G+A+T) < 14 , terminal_melting := 4*(C+G)+2*(A+T)]
    base_counts[(C+G+A+T) >= 14 , terminal_melting := 64.9 + 41*(G+C-16.4)/(A+T+G+C)]
 
    return(base_counts)
}

nearest_neighbors_terminal_melting_calculation <- function(sequence_list){
    require(rmelting)
    melting_list = foreach(seq = sequence_list) %do% {
        temp = melting(seq, nucleic.acid.conc = 2e-06, hybridisation.type = "dnadna", Na.conc=1) 
        unlist(temp)[names(unlist(temp)) == 'Results.Melting temperature (C)']
    }
    cleaned = unlist(melting_list)
    names(cleaned) = NULL
    column = data.table(terminal_melting = cleaned)
    return(column)
}

get_melting_temp <- function(calculation_type){
    stopifnot(calculation_type %in% c('simple', 'nearest_neighbors'))
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

    # get length of terminal seq
    #NOT INCLUDING PNUCS IN THIS TERMINAL GC CONTENT CALCULATION
    together[, depth := trim_length + LEFT_NUC_MOTIF_COUNT]

    # get terminal seq
    together[, terminal_seq := substring(sequences, nchar(sequences) - depth + 1, nchar(sequences))]

    if (calculation_type == 'simple'){
        seq_list = DNAStringSet(together$terminal_seq)
        melting_temps = simple_terminal_melting_calculation(seq_list) 
    } else if (calculation_type == 'nearest_neighbors'){
        seq_list = together$terminal_seq
        melting_temps =  nearest_neighbors_terminal_melting_calculation(seq_list) 
    }

    # merge
    together = cbind(together, melting_temps)
    return(unique(together[, c('gene', 'trim_length', 'terminal_melting')]))
}


process_for_terminal_melting <- function(group_motif_data, calculation_type){
    row_count = nrow(group_motif_data)
    if (!('terminal_melting' %in% colnames(group_motif_data))){
        terminal_gc = get_melting_temp(calculation_type)
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

