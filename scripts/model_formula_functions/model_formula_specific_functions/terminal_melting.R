simple_terminal_melting_calculation <- function(sequence_list){
    seq_list = DNAStringSet(sequence_list)
    # simple terminal melting calculation 
    base_counts = as.data.table(letterFrequency(seq_list, letters="ACGT", OR = 0))
    base_counts[(C+G+A+T) < 14 , terminal_melting := 4*(C+G)+2*(A+T)]
    base_counts[(C+G+A+T) >= 14 , terminal_melting := 64.9 + 41*(G+C-16.4)/(A+T+G+C)]
 
    temps = data.table(terminal_seq = sequence_list, terminal_melting = base_counts$terminal_melting)
    return(temps)
}

nearest_neighbors_terminal_melting_calculation <- function(sequence_list){
    require(rmelting)
    require(parallel)
    cluster = makeCluster(NCPU)
    melting_list = parLapply(cluster, sequence_list, function(x){
                                 melting = rmelting::melting(x, nucleic.acid.conc = 2e-06, hybridisation.type = "dnadna", Na.conc=1)
                                 temp = unlist(melting)[names(unlist(melting)) == 'Results.Melting temperature (C)']
                                 data.table::data.table(terminal_seq = x, terminal_melting = temp)
})
    stopCluster(cluster)
    
    cleaned = bind_rows(melting_list)
    cleaned$terminal_melting = as.numeric(cleaned$terminal_melting)
    return(cleaned)
}

combo_terminal_melting_calculation <- function(sequence_list){
    short_seqs = sequence_list[nchar(sequence_list) <= 8]
    long_seqs = sequence_list[nchar(sequence_list) > 8]
    
    # use simple melting temp calculation for short sequences and the nearest neighbors approach for long sequences
    short_temps = simple_terminal_melting_calculation(short_seqs)
    long_temps = nearest_neighbors_terminal_melting_calculation(long_seqs)

    return(rbind(short_temps, long_temps))
}

get_melting_temp <- function(calculation_type){
    stopifnot(calculation_type %in% c('simple', 'nearest_neighbors', 'combo'))
    whole_nucseq = get_oriented_whole_nucseqs()

    trims = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)
    
    genes = whole_nucseq$gene[substring(whole_nucseq$gene, 4, 4) == toupper(substring(GENE_NAME, 1, 1))]
    together = data.table(gene = rep(genes, length(trims)), trim_length = rep(trims, length(genes)))
    together = merge(together, whole_nucseq, by = 'gene')

    setnames(together, 'gene', GENE_NAME)
    map = get_common_genes_from_seqs(together)
    together = merge(together, map, by = GENE_NAME)
    cols = c('trim_length', 'sequence', 'gene')
    together = together[, ..cols]

    # get length of terminal seq
    together[, depth := trim_length + LEFT_NUC_MOTIF_COUNT]

    # get terminal seq (only interested in double stranded end (not including pnucs or pnuc pairs))
    together[, terminal_seq := substring(sequence, nchar(sequence) - depth + 1, nchar(sequence)-abs(PNUC_COUNT))]

    if (calculation_type == 'simple'){
        melting_temps = simple_terminal_melting_calculation(together$terminal_seq) 
    } else if (calculation_type == 'nearest_neighbors'){
        melting_temps =  nearest_neighbors_terminal_melting_calculation(together$terminal_seq) 
    } else if (calculation_type == 'combo'){
        melting_temps = combo_terminal_melting_calculation(together$terminal_seq)    
    }

    # merge
    together = merge(together, unique(melting_temps), by = 'terminal_seq')

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

