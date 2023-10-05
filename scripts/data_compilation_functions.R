source(paste0(MOD_PROJECT_PATH, '/scripts/motif_class_functions/', MOTIF_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/annotation_specific_functions/', ANNOTATION_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/gene_specific_functions/', TRIM_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_formula_functions/', MODEL_TYPE, '.R'))

REQUIRED_COMMON_NUCS_5 <<- UPPER_TRIM_BOUND + 10 

get_gene_order <- function(gene_type){
    genes = strsplit(gene_type, '_')[[1]][1]
    genes = strsplit(genes, '-')[[1]]
    return(c(paste0(genes, '_gene')))
}

get_trim_order <- function(trim_type){
    trims = strsplit(trim_type, '_')[[1]][1]
    trims = strsplit(trims, '-')[[1]]
    return(c(paste0(trims, '_trim')))
}

get_subject_motif_output_location <- function(){
    output_location = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'motif_data')
    return(output_location)
}

get_common_genes_from_seqs <- function(subject_data, gene_type = GENE_NAME){
    if ('gene' %in% colnames(subject_data)){
        setnames(subject_data, 'gene', gene_type)
    }
    gene_col = gene_type
    sequence_col = paste0(gene_type, '_sequence')
    cols = c(gene_col, sequence_col)
    subject_data = unique(subject_data[,..cols])
    terminal_length = UPPER_TRIM_BOUND + REQUIRED_COMMON_NUCS_5 
    subject_data[[paste0(gene_type, '_terminal_seq')]] = substring(subject_data[[sequence_col]], nchar(subject_data[[sequence_col]])-(terminal_length -1), nchar(subject_data[[sequence_col]]))
    subject_data[[paste0(gene_type, '_class')]] = str_split(subject_data[[gene_col]], fixed('*'), simplify = TRUE)[,1] 
    cols = c(paste0(gene_type, '_terminal_seq'), paste0(gene_type, '_class'))
    subject_data[,cdr3_gene_group := .GRP, by = cols]
    for (gene_class_group in unique(subject_data[[paste0(gene_type, '_class')]])){
        temp = subject_data[get(paste0(gene_type, '_class')) == gene_class_group]
        if (length(unique(temp$cdr3_gene_group)) == 1){
            subject_data[get(paste0(gene_type, '_class')) == gene_class_group, paste0('grouped_by_common_', gene_type, '_sequence') := TRUE]
        } else {
            subject_data[get(paste0(gene_type, '_class')) == gene_class_group, paste0(gene_type, '_class') := get(gene_col)]
            subject_data[get(paste0(gene_type, '_class')) == gene_class_group, paste0('grouped_by_common_', gene_type, '_sequence') := FALSE]
        }
    }
    return_cols = c(gene_col, paste0(gene_type, '_class'))
    data_subset = subject_data[, ..return_cols]
    setnames(data_subset, paste0(gene_type, '_class'), paste0(gene_type, '_group'))
    return(data_subset)
}

get_nuc_context <- function(whole_gene_seqs, trim_lengths){
     stopifnot(length(whole_gene_seqs) == length(trim_lengths))
     left_nuc_list = c()
     right_nuc_list = c()
     for (index in seq(1, length(whole_gene_seqs))){
        whole_gene_seq = DNAString(whole_gene_seqs[index])
        if (PNUC_COUNT > 0){
            possible_pnucs_5_to_3 = substring(reverseComplement(whole_gene_seq),1, PNUC_COUNT)
        } else if (PNUC_COUNT < 0){
            possible_pnucs_5_to_3 = DNAString() 
            whole_gene_seq = substring(whole_gene_seq, 1, nchar(whole_gene_seq) + PNUC_COUNT)             
        } else {
            possible_pnucs_5_to_3 = DNAString() 
        }

        whole_gene_and_pnucs = c(unlist(whole_gene_seq), unlist(possible_pnucs_5_to_3))

        trim_length = trim_lengths[index]
        trimmed_gene_seq = substring(whole_gene_and_pnucs, 1, nchar(whole_gene_seq)-trim_length)
        trimmed_length = nchar(trimmed_gene_seq)

        left_nucs = substring(trimmed_gene_seq, trimmed_length - (REQUIRED_COMMON_NUCS_5 - 1), trimmed_length) 

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

get_unobserved_nuc_context <- function(tcr_dataframe, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    # get observed trim lengths, genes
    remove_cols = c('count', paste0(genes, 'whole_seq'), paste0(trims, '_left_nucs'), paste0(trims, '_right_nucs'))
    cols = colnames(tcr_dataframe)[!(colnames(tcr_dataframe) %in% remove_cols)]
    tcr_dataframe_observed = tcr_dataframe[,.N, by = cols]
    
    # get unique genes
    non_trim_cols = cols[!(cols %in% c(trims, 'N'))]
    unique_obs = unique(tcr_dataframe_observed[,..non_trim_cols])

    # get desired trim lengths, genes
    trim_lengths = seq(LOWER_TRIM_BOUND,UPPER_TRIM_BOUND)

    desired_obs = unique_obs
    for (i in seq(length(trims))){
        desired_obs = desired_obs %>%
            mutate(!!trims[i]:= list(trim_lengths)) %>%
            unnest(cols = trims[i]) %>%
            as.data.table()
    }

    # get unobserved subset
    unobserved = desired_obs[!tcr_dataframe_observed, on = cols]

    # get unobserved nuc context
    for (i in seq(length(trims))){
        u_cols = c(paste0(genes[i], '_whole_seq'), paste0(genes[i], '_group'), trims[i])
        subset = unique(unobserved[, ..u_cols])
        subset = subset[, c(paste0(trims[i], '_left_nucs'), paste0(trims[i], '_right_nucs')):= get_nuc_context(get(paste0(genes[i], '_whole_seq')), get(trims[i]))] 
        unobserved = merge(unobserved, subset, by = u_cols)
    }

    unobserved$count = 0
    return(unobserved)
}

adaptive_data_filtering <- function(data){
    if ('adaptive_v_gene_call' %in% colnames(data)){
        data = data[v_gene == adaptive_v_gene_call]
        if (LOCUS == 'alpha'){
            data = data[vj_insert <= 15]
        }
    }
    return(data)
}

filter_by_productivity <- function(data){
    stopifnot(PRODUCTIVITY %in% c('productive', 'nonproductive', 'both'))
    if (PRODUCTIVITY == 'productive'){
        if (all(unique(data$productive) %in% c(TRUE, FALSE))){
            data = data[productive == TRUE]
        } else {
            data = data[productive == 'productive']
        }
    } else if (PRODUCTIVITY == 'nonproductive'){
        if (all(unique(data$productive) %in% c(TRUE, FALSE))){
            data = data[productive == FALSE]
        } else {
            data = data[productive == 'nonproductive']
        }
    }
    return(data)
}

condense_tcr_data <- function(tcr_dataframe, gene_type = GENE_NAME, trim_type = TRIM_TYPE, insertions = INSERTIONS){
    stopifnot(insertions %in% c('nonzero', 'zero', 'all'))
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    cols = c(paste0(genes, '_sequence'), paste0(genes, '_group'), paste(trims), JOINING_INSERT)
    if ('count' %in% colnames(tcr_dataframe)){
        cols = c(cols, 'count')
    }
    tcr_dataframe = tcr_dataframe[,..cols]
    if (insertions == 'nonzero'){
        tcr_dataframe = tcr_dataframe[get(JOINING_INSERT) != 0]
    } else if (insertions == 'zero'){
        tcr_dataframe = tcr_dataframe[get(JOINING_INSERT) == 0]
    }
    setnames(tcr_dataframe, paste0(genes, '_sequence'), paste0(genes, '_whole_seq'))
    #condense data by gene, trim, etc.
    cols = c(paste0(genes, '_whole_seq'), trims, paste0(genes, '_group'))
    if ('count' %in% colnames(tcr_dataframe)){
        condensed_tcr = tcr_dataframe[, sum(count), by = cols]
        setnames(condensed_tcr, 'V1', 'count')
    } else {
        condensed_tcr = tcr_dataframe[, .N, by = cols]
        setnames(condensed_tcr, 'N', 'count')
    }
    return(condensed_tcr)
}

sum_trim_observations <- function(condensed_tcr_dataframe, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    cols = c(paste0(genes, '_group'), trims, paste0(trims, '_left_nucs'), paste0(trims, '_right_nucs'))
    summed = condensed_tcr_dataframe[, sum(count), by = cols]
    setnames(summed, 'V1', 'count')
    return(summed)
}

general_get_all_nuc_contexts <- function(tcr_dataframe, subject_id, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    #filter data by trim bounds
    if (length(trims) > 1){
        tcr_dataframe = tcr_dataframe[get(trims[1]) >= LOWER_TRIM_BOUND & get(trims[1]) <= UPPER_TRIM_BOUND & get(trims[2]) >= LOWER_TRIM_BOUND & get(trims[2]) <= UPPER_TRIM_BOUND]
    } else {
         tcr_dataframe = tcr_dataframe[get(trims[1]) >= LOWER_TRIM_BOUND & get(trims[1]) <= UPPER_TRIM_BOUND]
    }
    if (nrow(tcr_dataframe) == 0){
        return(tcr_dataframe)
    } else {
        #condense data by gene, trim, etc.
        tcr_dataframe = condense_tcr_data(tcr_dataframe, gene_type = gene_type, trim_type = trim_type)

        #get motifs for observed trims
        motif_dataframe = tcr_dataframe
        for (i in seq(length(trims))){
            u_cols = c(paste0(genes[i], '_whole_seq'), paste0(genes[i], '_group'), trims[i])
            subset = unique(motif_dataframe[, ..u_cols])
            subset = subset[, c(paste0(trims[i], '_left_nucs'), paste0(trims[i], '_right_nucs')):= get_nuc_context(get(paste0(genes[i], '_whole_seq')), get(trims[i]))] 
            motif_dataframe = merge(motif_dataframe, subset, by = u_cols)
        }

        # get motifs for unobserved
        if (nrow(motif_dataframe) < length(interaction(motif_dataframe[[paste0(genes[1], '_group')]], motif_dataframe[[paste0(genes[2], '_group')]]))*(UPPER_TRIM_BOUND - LOWER_TRIM_BOUND + 1)*(UPPER_TRIM_BOUND - LOWER_TRIM_BOUND + 1)){
            unobserved = get_unobserved_nuc_context(motif_dataframe, gene_type = gene_type, trim_type = trim_type)
            together = rbind(motif_dataframe, unobserved)
        } else {
            together = motif_dataframe
        }

        # condense observations by group
        recondensed = sum_trim_observations(together, gene_type = gene_type, trim_type = trim_type)

        recondensed[count == 0, paste0(trim_type, '_observed') := FALSE]
        recondensed[count != 0, paste0(trim_type, '_observed') := TRUE]

        recondensed$gene_type = gene_type 
        if (!is.null(subject_id)){
            recondensed$subject = subject_id
        }
        return(recondensed)
    }
}

convert_data_to_motifs <- function(compiled_data, left_window_size = LEFT_NUC_MOTIF_COUNT, right_window_size = RIGHT_NUC_MOTIF_COUNT, trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)

    for (i in seq(length(trims))){
        left =paste0(trims[i], '_left_nucs') 
        right = paste0(trims[i], '_right_nucs')
        compiled_data[, paste0(trims[i], '_left_motif') := substring(get(left), REQUIRED_COMMON_NUCS_5-(left_window_size-1), REQUIRED_COMMON_NUCS_5)]
        compiled_data[, paste0(trims[i], '_right_motif') := substring(get(right), 1, right_window_size)]
        compiled_data[, paste0(trims[i], '_motif') := paste0(get(left), get(right))] 
    }
    return(compiled_data)
}

process_data_for_model_fit <- function(group_motif_data, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    processed = group_motif_data
    for (i in seq(length(trims))){
        processed = process_single_data_for_model_fit(processed, gene_type = genes[i], trim_type = trims[i])
    }
    return(processed)
}

split_motif_column_by_motif_position <- function(aggregated_subject_data, trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)

    for (i in seq(length(trims))){
        positions = get_positions(trims[i])
        for (side in c('5end', '3end')){
            subset = positions[positions %like% side]
            lr_side = ifelse(side == '5end', 'left', 'right')
            if (length(subset) > 0){
                aggregated_subject_data[, paste0(subset) := tstrsplit(get(paste0(trims[i], '_', lr_side, '_motif')), "")]
            } 
        }
    }
    return(aggregated_subject_data)
}

set_contrasts <- function(group_motif_data, ref_base = 'A', trim_type = TRIM_TYPE){
    trims = get_trim_order(trim_type)

    for (i in seq(length(trims))){
        group_motif_data = single_set_contrasts(group_motif_data, ref_base, trims[i])
    }
    return(group_motif_data)
}

compile_data_for_subject <- function(file_path, write = TRUE, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    temp_data = fread(file_path)
    temp_data = reformat_data(temp_data)
    if (gene_type == 'd_gene'){
        temp_data = temp_data[d_gene != '-']
    }

    subject_id = extract_subject_ID(file_path)

    temp_data = filter_by_productivity(temp_data)    
    temp_data = adaptive_data_filtering(temp_data)
    output_location = get_subject_motif_output_location() 
    dir.create(output_location, recursive = TRUE, showWarnings = FALSE)
    together = get_oriented_full_sequences(temp_data, gene_type = gene_type)
    if ('subject' %in% colnames(together)){
        motif_data = data.table()
        for (subject_i in unique(together$subject)){
            temp = get_all_nuc_contexts(together[subject == subject_i], subject_id = subject_i, gene_type = gene_type, trim_type = trim_type)
            motif_data = rbind(motif_data, temp)
        }
    } else {
        motif_data = get_all_nuc_contexts(together, subject_id, gene_type = gene_type, trim_type = trim_type)
    }

    motif_data = filter_motif_data_for_possible_sites(motif_data)

    if (isTRUE(write)){
        fwrite(motif_data, file = file.path(output_location, paste0(subject_id, '.tsv')), sep = '\t')
    } else {
        return(motif_data)
    }
}

compile_all_data <- function(directory, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    foreach(file = files) %dopar% {
        compile_data_for_subject(file, gene_type = gene_type, trim_type = trim_type)
        print(paste0(file))
    }
    stopImplicitCluster()
}

get_frames_for_pair <- function(v_seq, v_frame, v_trim, j_seq, j_frame, j_trim){
    # get processed V-gene sequence
    adjusted_v_seq = substring(v_seq, v_frame) 
    if (v_trim < 0){
        adjusted_v_seq = paste0(adjusted_v_seq, paste0(rep('X', -1*v_trim), collapse = ''), collapse = '')
    } else {
        adjusted_v_seq = substring(adjusted_v_seq, 1, nchar(adjusted_v_seq) - v_trim)
    }

    # get processed J-gene sequence
    j_extra = (nchar(j_seq) - (j_frame - 1))%%3
    adjusted_j_seq = substring(j_seq, 1, nchar(j_seq) - j_extra)
    if (j_trim < 0){
        adjusted_j_seq = paste0(paste0(rep('X', -1*j_trim), collapse = ''), adjusted_j_seq, collapse = '')
    } else {
        adjusted_j_seq = substring(adjusted_j_seq, j_trim + 1)
    }

    tog = paste0(adjusted_v_seq, adjusted_j_seq, collapse = '')
    frame = nchar(tog) %% 3
    if (frame == 0){
        frame_type = 'In'
    } else {
        frame_type = 'Out'
    }
    return(frame_type)
}

get_all_frames <- function(){
    frames = fread('https://raw.githubusercontent.com/phbradley/conga/master/conga/tcrdist/db/combo_xcr.tsv')[organism == 'human' & chain == substring(CHAIN_TYPE, 3, 3)]

    v = frames[region == 'V'][, c('id', 'frame', 'nucseq')]
    colnames(v) = c('v_gene', 'v_frame', 'v_seq')
    j = frames[region == 'J'][, c('id', 'frame', 'nucseq')]
    colnames(j) = c('j_gene', 'j_frame', 'j_seq')

    v$dummy = 1
    j$dummy = 1

    pairs = merge(v, j, by = 'dummy', allow.cartesian = TRUE)

    trims = data.table(v_trim = rep(seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND), each = (UPPER_TRIM_BOUND - LOWER_TRIM_BOUND) + 1), j_trim = rep(seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND), (UPPER_TRIM_BOUND - LOWER_TRIM_BOUND) + 1))
    trims$dummy = 1

    pairs = merge(pairs, trims, by = 'dummy', allow.cartesian = TRUE)[, -c('dummy')]

    pairs[, frame_type := apply(.SD, 1, function(x) get_frames_for_pair(x[1], as.numeric(x[2]), as.numeric(x[3]), x[4], as.numeric(x[5]), as.numeric(x[6]))), .SDcols = c("v_seq", "v_frame", "v_trim", "j_seq", "j_frame", "j_trim")]

    return(pairs[, -c('v_seq', 'j_seq')])
}

get_frames_data <- function(){
    path = file.path(ROOT_PATH, 'data')
    dir.create(path)
    file_name = file.path(path, paste0(CHAIN_TYPE, '_paired_frame_data_trims_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND, '.tsv'))

    if (!file.exists(file_name)) {
        frame_data = get_all_frames()
        fwrite(frame_data, file_name, sep = '\t')
    } else {
        frame_data = fread(file_name)
    }
    return(frame_data)
}

processed_data_path <- function(){
    output_location = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), paste0(LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_', MODEL_TYPE))
    dir.create(output_location, recursive = TRUE, showWarnings = FALSE)
    filename = file.path(output_location, 'processed_data.tsv')
    return(filename)
}

subset_processed_data <- function(data, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    params = get_parameter_vector(trims, genes)
    other = c(paste0(genes, '_group'), trims, 'weighted_observation', 'count', 'total_tcr')
    cols = c(other, params)
    return(data[, ..cols])
}
