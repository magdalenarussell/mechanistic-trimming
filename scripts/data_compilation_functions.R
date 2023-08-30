source(paste0(MOD_PROJECT_PATH, '/scripts/motif_class_functions/', MOTIF_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/annotation_specific_functions/', ANNOTATION_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/gene_specific_functions/', TRIM_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_formula_functions/', MODEL_TYPE, '.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/data_grouping_functions/', DATA_GROUP, '.R'))

REQUIRED_COMMON_NUCS_5 <<- 10

get_subject_motif_output_location <- function(){
    output_location = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, DATA_GROUP, TRIM_TYPE, PRODUCTIVITY, paste0(MOTIF_TYPE, '_motif_trims_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND), 'motif_data')
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

get_oriented_full_sequences <- function(subject_data, whole_nucseq = get_oriented_whole_nucseqs(), gene_type = GENE_NAME){
    temp_data = merge(subject_data, whole_nucseq, by.x = gene_type, by.y = 'gene')
    gene_seqs = whole_nucseq[substring(gene, 4, 4) == toupper(substring(gene_type, 1, 1))]
    gene_groups = get_common_genes_from_seqs(gene_seqs, gene_type)
    together = merge(temp_data, gene_groups, by = gene_type)
    return(together)
}

get_unobserved_nuc_context <- function(tcr_dataframe, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    # get observed trim lengths, genes
    remove_cols = c('count', paste0(trim_type, 'med_seq'), paste0(gene_type, 'whole_seq'), paste0(trim_type, '_left_nucs'), paste0(trim_type, '_right_nucs'), paste0(trim_type, '_observed'))
    cols = colnames(tcr_dataframe)[!(colnames(tcr_dataframe) %in% remove_cols)]
    tcr_dataframe_observed = tcr_dataframe[,.N, by = cols]
    # get unique genes
    non_trim_cols = cols[!(cols %in% c(trim_type, 'N'))]
    unique_obs = unique(tcr_dataframe_observed[,..non_trim_cols])
    # get desired trim lengths, genes
    trim_lengths = seq(LOWER_TRIM_BOUND,UPPER_TRIM_BOUND)
    desired_obs = unique_obs %>%
        mutate(!!trim_type := list(trim_lengths)) %>%
        unnest(cols = c(trim_type)) %>%
        as.data.table()
    # get unobserved subset
    unobserved = desired_obs[!tcr_dataframe_observed, on = cols]
    # get unobserved nuc context
    unobserved = unobserved[, c(paste0(trim_type, '_left_nucs'), paste0(trim_type, '_right_nucs')):= get_nuc_context(get(paste0(gene_type, '_whole_seq')), get(trim_type))] 
    unobserved[[paste0(trim_type, '_observed')]] = FALSE
    unobserved$count = 0
    cols = colnames(unobserved)[!(colnames(unobserved) %like% 'whole_seq')]
    return(unobserved[,..cols])
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

general_get_all_nuc_contexts <- function(tcr_dataframe, subject_id, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    #filter data by trim bounds
    tcr_dataframe = tcr_dataframe[get(TRIM_TYPE) >= LOWER_TRIM_BOUND & get(TRIM_TYPE) <= UPPER_TRIM_BOUND]
    if (nrow(tcr_dataframe) == 0){
        return(tcr_dataframe)
    } else {
        #condense data by gene, trim, etc.
        tcr_dataframe[, paste0(trim_type, 'med_seq') := substring(get(paste0(gene_type, '_sequence')), 1, nchar(get(paste0(gene_type, '_sequence'))) - get(trim_type))]
        tcr_dataframe = condense_tcr_data(tcr_dataframe, gene_type = gene_type, trim_type = trim_type)
        #get motifs
        motif_dataframe = tcr_dataframe[,c(paste0(trim_type, '_left_nucs'), paste0(trim_type, '_right_nucs')):=get_nuc_context(get(paste0(gene_type, '_whole_seq')), get(trim_type))]
        motif_dataframe[[paste0(trim_type, '_observed')]] = TRUE
        if (nrow(motif_dataframe) < length(unique(motif_dataframe[[paste0(gene_type, '_group')]]))*(UPPER_TRIM_BOUND - LOWER_TRIM_BOUND + 1)){
            unobserved = get_unobserved_nuc_context(motif_dataframe, gene_type = gene_type, trim_type = trim_type)
            together = rbind(motif_dataframe[,-c(1,2)], unobserved)
        } else {
            together = motif_dataframe[,-c(1,2)]
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

get_colnames <- function(){
    names = list(v_call="v_call", d_call="d_call", j_call="j_call", v_germline_start="v_germ_start_vdj", d_germline_start="d_germ_start", j_germline_start="j_germ_start", np1_length="np1_length", np2_length="np2_length", junction="junction", junction_length="junction_length", sequence_alignment="sequence_alignment")
    return(names)
}

count_deletions_changeo <- function(row, type, gene_type = GENE_NAME){
    stopifnot(type %in% c('v', 'j'))
    gene_type = paste0(substring(type, 1, 1), '_call')
    allele = str_split(row[[gene_type]], ',')[[1]][1]
    deleted = c(NA, NA)
    if (is.na(allele)) { 
        return(deleted) 
    }

    germline_seqs = get_whole_nucseqs() 
    germline = germline_seqs[gene == allele]$sequence

    names = get_colnames()
    allele_germline_start = as.numeric(row[[as.character(names[paste0(type, '_germline_start')])]])
    allele_germline_end = allele_germline_start + row[[paste0(type, '_seq_length')]] - 1
        
    germline_head = stringi::stri_sub(germline, 1, allele_germline_start - 1)
    deleted_head = nchar(gsub("\\.", "", germline_head))
    germline_tail = stringi::stri_sub(germline, allele_germline_end+1, nchar(germline))
    deleted_tail = nchar(gsub("\\.", "", germline_tail))
                
    deleted[1] = deleted_head
    deleted[2] = deleted_tail
    return(deleted)
}

sample_by_family <- function(data){
    counts = data[, .N, by = family]
    singles = counts[N == 1]$family
    not_single = counts[N != 1]$family
    single_data = data[family %in% singles]
    for (fam in not_single){
        max = data[family == fam, max(conscount)]
        selected = data[family == fam & conscount == max]
        single_data = rbind(single_data, selected)
    }
    return(single_data)
}

process_changeo <- function(file_path, gene_type = GENE_NAME){
    file_data = fread(file_path)
    colnames(file_data) = tolower(colnames(file_data))
    for (i in 1:nrow(file_data)){
        row = file_data[i]
        v_dels = count_deletions_changeo(row, type = 'v', gene_type = gene_type) 
        j_dels = count_deletions_changeo(row, type = 'j', gene_type = gene_type)
        file_data[['v_trim']][i] = v_dels[2]
        file_data[['j_trim']][i] = j_dels[1]
        file_data[['v_gene']][i] = str_split(row[['v_call']], ',')[[1]][1]
        file_data[['j_gene']][i] = str_split(row[['j_call']], ',')[[1]][1]
    }
    file_data = sample_by_family(file_data)
    cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim')
    subset = file_data[, ..cols]
    subset$productive = FALSE
    return(subset)
}

compile_data_for_subject <- function(file_path, write = TRUE, gene_type = GENE_NAME, trim_type = TRIM_TYPE){
    if (file_path %like% 'tab'){
        temp_data = process_changeo(file_path, gene_type)
    } else {
        temp_data = fread(file_path)
    }
    if (gene_type == 'd_gene'){
        temp_data = temp_data[d_gene != '-']
    }
    if (!('subject' %in% colnames(temp_data))){
        subject_id = extract_subject_ID(file_path)
    }else {
        subject_id = NULL
    }
    temp_data = filter_by_productivity(temp_data)    
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

convert_data_to_motifs <- function(compiled_data, left_window_size = LEFT_NUC_MOTIF_COUNT, right_window_size = RIGHT_NUC_MOTIF_COUNT, trim_type = TRIM_TYPE){
    left =paste0(trim_type, '_left_nucs') 
    right = paste0(trim_type, '_right_nucs')
    compiled_data[, paste0(trim_type, '_left_motif') := substring(get(left), REQUIRED_COMMON_NUCS_5-(left_window_size-1), REQUIRED_COMMON_NUCS_5)]
    compiled_data[, paste0(trim_type, '_right_motif') := substring(get(right), 1, right_window_size)]
    compiled_data[, paste0(trim_type, '_motif') := paste0(get(left), get(right))] 
    return(compiled_data)
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
