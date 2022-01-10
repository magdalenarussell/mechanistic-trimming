get_all_base_variables <- function(side){
    stopifnot(side %in% c('left', 'right'))
    bases = c('A', 'T', 'G', 'C')
    vars = paste0(side, '_base_count_', bases)
    return(vars)
}

count_bases_seq_list <- function(seq_list, side){
    subset = unique(seq_list)
    bases = c('A', 'T', 'G', 'C')
    registerDoParallel(cores=NCPU)  
    cl = makeCluster(NCPU, type="FORK")
    seq_list_DNA = parLapply(cl, subset, function(x) DNAString(x))
    counts = parLapply(cl, seq_list_DNA, function(x) letterFrequency(x, letters = bases))
    counts_dt = rbindlist(lapply(counts, as.data.frame.list))
    stopCluster(cl)  
    for (base in bases) {
        setnames(counts_dt, base, paste0(side, '_base_count_', base))
    }
    together = cbind(subset, counts_dt)
    setnames(together, 'subset', paste0(side, '_motif_nucs'))
    return(together)
}
    
process_for_two_side_base_count <- function(motif_data){
    vars = get_all_base_variables('left')
    if (!all(vars %in% colnames(motif_data))){
        motif_data[, left_motif_nucs := substring(motif, 1, LEFT_NUC_MOTIF_COUNT)]
        motif_data[, right_motif_nucs := substring(motif, LEFT_NUC_MOTIF_COUNT + 1, RIGHT_NUC_MOTIF_COUNT +LEFT_NUC_MOTIF_COUNT)]

        for (side in c('left', 'right')){
            col = paste0(side, '_motif_nucs')
            base_counts = count_bases_seq_list(motif_data[[col]], side)
            motif_data = merge(motif_data, base_counts, by = col)
        }
    }
    return(motif_data)
}
