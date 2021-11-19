source(paste0('scripts/motif_class_functions/', MOTIF_TYPE, '.R'))

extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][3]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_")[[1]][3]
    return(localID)
}

get_subject_motif_output_location <- function(){
    output_location = file.path(OUTPUT_PATH, TRIM_TYPE, paste0(MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
    return(output_location)
}

get_common_genes_from_seqs <- function(subject_data){
    cols = c(GENE_NAME, 'sequences')
    subject_data = unique(subject_data[,..cols])
    terminal_length = UPPER_TRIM_BOUND + LEFT_NUC_MOTIF_COUNT
    subject_data$terminal_seq = substring(subject_data$sequences, nchar(subject_data$sequences)-(terminal_length -1), nchar(subject_data$sequences))
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

#TODO add 5' or 3' designation...for orientation
get_motif_context <- function(whole_gene_seqs, trimmed_gene_seqs, trim_lengths){
    stopifnot(length(whole_gene_seqs) == length(trimmed_gene_seqs))
    stopifnot(length(whole_gene_seqs) == length(trim_lengths))
    motifs = c()
    for (index in seq(1, length(whole_gene_seqs))){
        whole_gene_seq = DNAString(whole_gene_seqs[index])
        trimmed_gene_seq = DNAString(trimmed_gene_seqs[index])
        trimmed_length = nchar(trimmed_gene_seqs[index])
        original_trimmed_length = trimmed_length
        trim_length = trim_lengths[index]

        if (nchar(trimmed_gene_seq) < LEFT_NUC_MOTIF_COUNT){
            trimmed_gene_seq = substring(whole_gene_seq, 1, nchar(whole_gene_seq)-trim_length)
            trimmed_length = nchar(trimmed_gene_seq)
        }

        left_nuc_motif = substring(trimmed_gene_seq, trimmed_length - (LEFT_NUC_MOTIF_COUNT - 1), trimmed_length) 
        
        if (PNUC_COUNT > 0){
            possible_pnucs_5_to_3 = substring(reverseComplement(whole_gene_seq),1, PNUC_COUNT)
        } else {
            possible_pnucs_5_to_3 = DNAString() 
        }
        whole_gene_and_pnucs = c(unlist(whole_gene_seq), unlist(possible_pnucs_5_to_3))
        seq_right_of_trim = substring(whole_gene_and_pnucs, nchar(whole_gene_and_pnucs)-(trim_length + PNUC_COUNT) + 1, nchar(whole_gene_and_pnucs))

        if (nchar(seq_right_of_trim) < RIGHT_NUC_MOTIF_COUNT){
            missing_nucs = DNAString(strrep('-', RIGHT_NUC_MOTIF_COUNT - nchar(seq_right_of_trim)))
            seq_right_of_trim = c(unlist(seq_right_of_trim), unlist(missing_nucs))
        }

        right_nuc_motif = substring(seq_right_of_trim, 1, RIGHT_NUC_MOTIF_COUNT)

        motif = c(unlist(left_nuc_motif), unlist(right_nuc_motif))
        motifs = c(motifs, as.character(motif))
    }
    return(motifs)
}

#TODO add 5' or 3' designation...for orientation
get_motif_context_unobserved <- function(whole_gene_seqs, trim_lengths){
    stopifnot(length(whole_gene_seqs) == length(trim_lengths))
    motifs = c()
    for (index in seq(1, length(whole_gene_seqs))){
        whole_gene_seq = DNAString(whole_gene_seqs[index])
        trim_length = trim_lengths[index]

        trimmed_gene_seq = substring(whole_gene_seq, 1, nchar(whole_gene_seq)-trim_length)
        trimmed_length = nchar(trimmed_gene_seq)

        left_nuc_motif = substring(trimmed_gene_seq, trimmed_length - (LEFT_NUC_MOTIF_COUNT - 1), trimmed_length) 

        if (PNUC_COUNT > 0){
            possible_pnucs_5_to_3 = substring(reverseComplement(whole_gene_seq),1, PNUC_COUNT)
        } else {
            possible_pnucs_5_to_3 = DNAString() 
        }
        whole_gene_and_pnucs = c(unlist(whole_gene_seq), unlist(possible_pnucs_5_to_3))
        seq_right_of_trim = substring(whole_gene_and_pnucs, nchar(whole_gene_and_pnucs)-(trim_length + PNUC_COUNT) + 1, nchar(whole_gene_and_pnucs))

        if (nchar(seq_right_of_trim) < RIGHT_NUC_MOTIF_COUNT){
            missing_nucs = DNAString(strrep('-', RIGHT_NUC_MOTIF_COUNT - nchar(seq_right_of_trim)))
            seq_right_of_trim = c(unlist(seq_right_of_trim), unlist(missing_nucs))
        }

        right_nuc_motif = substring(seq_right_of_trim, 1, RIGHT_NUC_MOTIF_COUNT)

        motif = c(unlist(left_nuc_motif), unlist(right_nuc_motif))
        motifs = c(motifs, as.character(motif))
    }
    return(motifs)
}

get_unobserved_motifs <- function(tcr_dataframe){
    tcr_dataframe_observed = tcr_dataframe[,.N, by = .(gene, trim_length)]
    unobserved_df = data.frame()
    for (gene_name in unique(tcr_dataframe_observed$gene)){
        observed = unique(tcr_dataframe_observed[gene == paste(gene_name)]$trim_length)
        unobserved = setdiff(seq(LOWER_TRIM_BOUND,UPPER_TRIM_BOUND), observed)
        unobserved_df = rbind(unobserved_df, data.frame(gene = rep(gene_name, length(unobserved)), trim_length = unobserved))
    }
    cols = c('whole_seq', 'gene', GENE_NAME)
    tcr_dataframe = unique(tcr_dataframe[,..cols])
    together = as.data.table(merge(unobserved_df, tcr_dataframe, by = 'gene'))
    together = together[, motif:= get_motif_context_unobserved(whole_seq, trim_length)] 
    together$observed = FALSE
    together$count = 0
    return(together[,-c('whole_seq')])
}

general_get_motifs <- function(tcr_dataframe, subject_id){
    #filter data by trim bounds
    tcr_dataframe = tcr_dataframe[get(TRIM_TYPE) >= LOWER_TRIM_BOUND & get(TRIM_TYPE) <= UPPER_TRIM_BOUND]
    cdr3_variable = paste0('cdr3_nucseq_from_', substring(TRIM_TYPE, 1, 1))
    cols = c(paste(cdr3_variable), 'sequences', 'gene', paste(TRIM_TYPE), paste0(GENE_NAME))
    tcr_dataframe = tcr_dataframe[,..cols]
    colnames(tcr_dataframe) = c('trimmed_seq', 'whole_seq', 'gene', 'trim_length', GENE_NAME)
    #condense data by gene, trim, etc.
    condensed_tcr = tcr_dataframe[, .N, by = .(trimmed_seq, whole_seq, gene, trim_length, get(GENE_NAME))]
    colnames(condensed_tcr) = c('trimmed_seq', 'whole_seq', 'gene', 'trim_length', GENE_NAME, 'count')
    #get motifs
    motif_dataframe = condensed_tcr[,motif:=get_motif_context(whole_seq, trimmed_seq, trim_length)]
    motif_dataframe$observed = TRUE
    unobserved = get_unobserved_motifs(motif_dataframe)
    together = rbind(motif_dataframe[,-c(1,2)], unobserved)
    together$gene_type = GENE_NAME
    together$subject = subject_id
    return(together)
}

compile_motifs_for_subject <- function(file_path){
    temp_data = fread(file_path)
    if (GENE_NAME == 'd_gene'){
        temp_data = temp_data[d_gene != '-']
    }
    subject_id = extract_subject_ID(file_path)
    
    output_location = get_subject_motif_output_location() 
    dir.create(output_location, recursive = TRUE, showWarnings = FALSE)
    whole_nucseq = fread('_ignore/tcrb_processed_geneseq.tsv')
    temp_data = merge(temp_data, whole_nucseq, by.x = GENE_NAME, by.y = 'gene_names')
    gene_seqs = whole_nucseq[substring(gene_names, 4, 4) == toupper(substring(GENE_NAME, 1, 1))]
    colnames(gene_seqs) = c(GENE_NAME, 'sequences')
    gene_groups = get_common_genes_from_seqs(gene_seqs)
    together = merge(temp_data, gene_groups, by = GENE_NAME)

    motif_data = get_motifs(together, subject_id)

    fwrite(motif_data, file = file.path(output_location, paste0(subject_id, '.tsv')), sep = '\t')
}

compile_all_motifs <- function(directory){
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    foreach(file = files) %dopar% {
        compile_motifs_for_subject(file)
        print(paste0(file))
    }
    stopImplicitCluster()
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
