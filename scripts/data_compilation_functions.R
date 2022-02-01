source(paste0('scripts/motif_class_functions/', MOTIF_TYPE, '.R'))
source(paste0('scripts/annotation_specific_functions/', ANNOTATION_TYPE, '.R'))
source(paste0('scripts/gene_specific_functions/', TRIM_TYPE, '.R'))
source(paste0('scripts/model_formula_functions/', MODEL_TYPE, '.R'))
source(paste0('scripts/data_grouping_functions/', DATA_GROUP, '.R'))

REQUIRED_COMMON_NUCS_5 <<- 10

get_subject_motif_output_location <- function(){
    output_location = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'motif_data', DATA_GROUP)
    return(output_location)
}

get_common_genes_from_seqs <- function(subject_data){
    cols = c(GENE_NAME, 'sequence')
    subject_data = unique(subject_data[,..cols])
    terminal_length = UPPER_TRIM_BOUND + REQUIRED_COMMON_NUCS_5 
    subject_data$terminal_seq = substring(subject_data$sequence, nchar(subject_data$sequence)-(terminal_length -1), nchar(subject_data$sequence))
    subject_data$gene_class = str_split(subject_data[[GENE_NAME]], fixed('*'), simplify = TRUE)[,1] 
    subject_data[,cdr3_gene_group := .GRP, by = .(terminal_seq, gene_class)]
    for (gene_class_group in unique(subject_data$gene_class)){
        temp = subject_data[gene_class == gene_class_group]
        if (length(unique(temp$cdr3_gene_group)) == 1){
            subject_data[gene_class == gene_class_group, grouped_by_common_sequence := TRUE]
        } else {
            subject_data[gene_class == gene_class_group, gene_class := get(GENE_NAME)]
            subject_data[gene_class == gene_class_group, grouped_by_common_sequence := FALSE]
        }
    }
    return_cols = c(GENE_NAME, 'gene_class')
    data_subset = subject_data[, ..return_cols]
    setnames(data_subset, 'gene_class', 'gene')
    return(data_subset)
}

get_nuc_context <- function(whole_gene_seqs, trim_lengths){
     stopifnot(length(whole_gene_seqs) == length(trim_lengths))
     left_nuc_list = c()
     right_nuc_list = c()
     for (index in seq(1, length(whole_gene_seqs))){
        whole_gene_seq = DNAString(whole_gene_seqs[index])
        trim_length = trim_lengths[index]
        trimmed_gene_seq = substring(whole_gene_seq, 1, nchar(whole_gene_seq)-trim_length)
        trimmed_length = nchar(trimmed_gene_seq)
        original_trimmed_length = trimmed_length

        left_nucs = substring(trimmed_gene_seq, trimmed_length - (REQUIRED_COMMON_NUCS_5 - 1), trimmed_length) 
        if (PNUC_COUNT > 0){
            possible_pnucs_5_to_3 = substring(reverseComplement(whole_gene_seq),1, PNUC_COUNT)
        } else if (PNUC_COUNT < 0){
            possible_pnucs_5_to_3 = DNAString() 
            whole_gene_seq = substring(whole_gene_seq, 1, nchar(whole_gene_seq) + PNUC_COUNT)             
        } else {
            possible_pnucs_5_to_3 = DNAString() 
        }

        whole_gene_and_pnucs = c(unlist(whole_gene_seq), unlist(possible_pnucs_5_to_3))
        seq_right_of_trim = substring(whole_gene_and_pnucs, nchar(whole_gene_and_pnucs)-(trim_length + PNUC_COUNT) + 1, nchar(whole_gene_and_pnucs))

        if (nchar(seq_right_of_trim) < REQUIRED_COMMON_NUCS_5){
            missing_nucs = DNAString(strrep('-', REQUIRED_COMMON_NUCS_5 - nchar(seq_right_of_trim)))
            seq_right_of_trim = c(unlist(seq_right_of_trim), unlist(missing_nucs))
        } else if (nchar(seq_right_of_trim) > REQUIRED_COMMON_NUCS_5){
            seq_right_of_trim = substring(seq_right_of_trim, 1, REQUIRED_COMMON_NUCS_5)
        }

        left_nuc_list = c(left_nuc_list, as.character(left_nucs))
        right_nuc_list = c(right_nuc_list, as.character(seq_right_of_trim))
    }
    return(list(left_nucs = left_nuc_list, right_nucs = right_nuc_list))
}

get_oriented_full_sequences <- function(subject_data){
    whole_nucseq = get_oriented_whole_nucseqs()
    temp_data = merge(subject_data, whole_nucseq, by.x = GENE_NAME, by.y = 'gene')
    gene_seqs = whole_nucseq[substring(gene, 4, 4) == toupper(substring(GENE_NAME, 1, 1))]
    setnames(gene_seqs, 'gene', GENE_NAME)
    gene_groups = get_common_genes_from_seqs(gene_seqs)
    together = merge(temp_data, gene_groups, by = GENE_NAME)
    return(together)
}

get_unobserved_nuc_context <- function(tcr_dataframe){
    # get observed trim lengths, genes
    remove_cols = c('count', 'trimmed_seq', 'whole_seq', 'left_nucs', 'right_nucs', 'observed')
    cols = colnames(tcr_dataframe)[!(colnames(tcr_dataframe) %in% remove_cols)]
    tcr_dataframe_observed = tcr_dataframe[,.N, by = cols]
    # get unique genes
    non_trim_cols = cols[!(cols %in% c('trim_length', 'N'))]
    unique_obs = unique(tcr_dataframe_observed[,..non_trim_cols])
    # get desired trim lengths, genes
    trim_lengths = seq(LOWER_TRIM_BOUND,UPPER_TRIM_BOUND)
    desired_obs = unique_obs %>%
        mutate(trim_length = list(trim_lengths)) %>%
        unnest(cols = c(trim_length)) %>%
        as.data.table()
    # get unobserved subset
    unobserved = desired_obs[!tcr_dataframe_observed, on = cols]
    # get unobserved nuc context
    new_cols = c('whole_seq', 'gene', GENE_NAME)
    tcr_dataframe = unique(tcr_dataframe[,..new_cols])
    together = as.data.table(merge(unobserved, tcr_dataframe, by = c(GENE_NAME, 'gene')))
    together = together[, c('left_nucs', 'right_nucs'):= get_nuc_context(whole_seq, trim_length)] 
    together$observed = FALSE
    together$count = 0
    return(together[,-c('whole_seq')])
}

filter_by_productivity <- function(data){
    stopifnot(PRODUCTIVITY %in% c('productive', 'nonproductive', 'both'))
    if (PRODUCTIVITY == 'productive'){
        data = data[productive == TRUE]
    } else if (PRODUCTIVITY == 'nonproductive'){
        data = data[productive == FALSE]
    }
    return(data)
}

general_get_all_nuc_contexts <- function(tcr_dataframe, subject_id){
    #filter data by trim bounds
    tcr_dataframe = tcr_dataframe[get(TRIM_TYPE) >= LOWER_TRIM_BOUND & get(TRIM_TYPE) <= UPPER_TRIM_BOUND]
    #condense data by gene, trim, etc.
    tcr_dataframe[, trimmed_seq := substring(sequence, 1, nchar(sequence) - get(TRIM_TYPE))]
    tcr_dataframe = condense_tcr_data(tcr_dataframe)
    #get motifs
    motif_dataframe = tcr_dataframe[,c('left_nucs', 'right_nucs'):=get_nuc_context(whole_seq, trim_length)]
    motif_dataframe$observed = TRUE
    unobserved = get_unobserved_nuc_context(motif_dataframe)
    together = rbind(motif_dataframe[,-c(1,2)], unobserved)
    # condense observations by group
    recondensed = sum_trim_observations(together)
    recondensed[count == 0, observed := FALSE]
    recondensed[count != 0, observed := TRUE]

    recondensed$gene_type = GENE_NAME
    recondensed$subject = subject_id
    return(recondensed)
}

compile_data_for_subject <- function(file_path){
    temp_data = fread(file_path)
    if (GENE_NAME == 'd_gene'){
        temp_data = temp_data[d_gene != '-']
    }
    subject_id = extract_subject_ID(file_path)
    temp_data = filter_by_productivity(temp_data)    
    output_location = get_subject_motif_output_location() 
    dir.create(output_location, recursive = TRUE, showWarnings = FALSE)
    together = get_oriented_full_sequences(temp_data)
    motif_data = get_all_nuc_contexts(together, subject_id)

    fwrite(motif_data, file = file.path(output_location, paste0(subject_id, '.tsv')), sep = '\t')
}

compile_all_data <- function(directory){
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    foreach(file = files) %dopar% {
        compile_data_for_subject(file)
        print(paste0(file))
    }
    stopImplicitCluster()
}

convert_data_to_motifs <- function(compiled_data){
    compiled_data[, left_motif := substring(left_nucs, REQUIRED_COMMON_NUCS_5-(LEFT_NUC_MOTIF_COUNT-1), REQUIRED_COMMON_NUCS_5)]
    compiled_data[, right_motif := substring(right_nucs, 1, RIGHT_NUC_MOTIF_COUNT)]
    compiled_data[, motif := paste0(left_motif, right_motif)] 
    return(compiled_data[, -c('left_nucs', 'right_nucs', 'left_motif', 'right_motif')])
}

get_background_freq_by_postion <- function(motif_data){
    positions = get_positions()
    backgrounds = data.table()
    for (pos in positions){
        back = motif_data[, sum(weighted_observation), by =pos]
        colnames(back) = c('base', 'count')
        back$total = sum(back$count)
        back$position = pos
        back$freq = back$count / back$total
        backgrounds = rbind(backgrounds, back)
    }
    return(backgrounds)
}

get_raw_cdr3_seqs <- function(tcr_repertoire_file){
    data = fread(tcr_repertoire_file)
    cdr3s = data$cdr3_nucseq
    return(cdr3s)
}
