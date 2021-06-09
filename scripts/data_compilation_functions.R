source(paste0('scripts/motif_class_functions/', MOTIF_TYPE, '.R'))

extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][3]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_")[[1]][3]
    return(localID)
}

get_subject_motif_output_location <- function(){
    output_location = file.path(OUTPUT_PATH, TRIM_TYPE, paste0(MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT))
    return(output_location)
}

get_common_genes_from_seqs <- function(subject_data){
    cols = c(GENE_NAME, 'sequences')
    subject_data = unique(subject_data[,..cols])
    subject_data$first_22_seq = substring(subject_data$sequences, nchar(subject_data$sequences)-21, nchar(subject_data$sequences))
    subject_data$gene_class = str_split(subject_data[[GENE_NAME]], fixed('*'), simplify = TRUE)[,1] 
    subject_data[,cdr3_gene_group := .GRP, by = .(first_22_seq, gene_class)]
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

        possible_pnucs_5_to_3 = substring(reverseComplement(whole_gene_seq),1, 2)
        whole_gene_and_pnucs = c(unlist(whole_gene_seq), unlist(possible_pnucs_5_to_3))
        seq_right_of_trim = substring(whole_gene_and_pnucs, nchar(whole_gene_and_pnucs)-(trim_length + 2) + 1, nchar(whole_gene_and_pnucs))

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

        possible_pnucs_5_to_3 = substring(reverseComplement(whole_gene_seq),1, 2)
        whole_gene_and_pnucs = c(unlist(whole_gene_seq), unlist(possible_pnucs_5_to_3))
        seq_right_of_trim = substring(whole_gene_and_pnucs, nchar(whole_gene_and_pnucs)-(trim_length + 2) + 1, nchar(whole_gene_and_pnucs))

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
    gene_groups = get_common_genes_from_seqs(temp_data)
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


