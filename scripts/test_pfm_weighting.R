library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)

TRIM_TYPE <<- 'v_trim' 
stopifnot(TRIM_TYPE == 'v_trim')

MOTIF_TYPE <<- 'pnuc_motif' 
stopifnot(MOTIF_TYPE %in% c('pnuc_motif', 'no_pnuc_motif'))

PFM_TYPE <<- 'unbounded'
stopifnot(PFM_TYPE %in% c('unbounded'))

NCPU <<- 1

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 4
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 4

source('scripts/pfm_functions.R')

get_gene_usage_frequency_by_subject <- function(gene, subject_id){
    usage = fread('/fh/fast/matsen_e/shared/tcr-gwas/exomotif/empirical_gene_frequencies.tsv')
    return(usage[subject == subject_id][[gene]])
}

get_gene_from_file_name <- function(file){
    file_name = strsplit(file, '/')[[1]][5]
    subset = strsplit(file_name, '_')[[1]][3]
    gene = strsplit(subset, '.tsv')[[1]][1]
    return(gene)
}

compile_PFMs_old <- function(group_localID_list, group_name, weighting = NULL){
    compiled_PFM = matrix(0, ncol = (RIGHT_NUC_MOTIF_COUNT+LEFT_NUC_MOTIF_COUNT), nrow = 5)
    for (subject in group_localID_list){
        path = file.path('_ignore', 'pfms', paste0(PFM_TYPE, '_', MOTIF_TYPE), paste0(subject))
        files = fs::dir_ls(path = path)
        for (file in files){
            if (!is.null(weighting)){
                #TODO add weighting scheme
            }
            temp_matrix = as.matrix(fread(file), rownames = 1)
            compiled_PFM = compiled_PFM + temp_matrix
        }
    }
    return(compiled_PFM)
}

compile_PPMs_new <- function(group_localID_list, group_name, weighting = NULL){
    compiled_PFM = matrix(0, ncol = (RIGHT_NUC_MOTIF_COUNT+LEFT_NUC_MOTIF_COUNT), nrow = 5)
    for (subject in group_localID_list){
        path = file.path('_ignore', 'pfms', paste0(PFM_TYPE, '_', MOTIF_TYPE), paste0(subject))
        files = fs::dir_ls(path = path)
        for (file in files){
            if (!is.null(weighting)){
                #TODO add weighting scheme
            }
            gene = get_gene_from_file_name(file) 
            temp_matrix = as.matrix(fread(file), rownames = 1)
            if (temp_matrix #TODO)
            subject_gene_usage = get_gene_usage_frequency_by_subject(gene, subject)
            gene_PFM = t(t(temp_matrix)/colSums(temp_matrix))
            gene_PFM_weighted = gene_PFM * subject_gene_usage
            compiled_PFM = compiled_PFM + gene_PFM_weighted
        }
    }
    return(compiled_PFM)
}

find_gene_nt_frequencies <- function(gene){
    whole_nucseq = fread('_ignore/tcrb_processed_geneseq.tsv')
    whole_nucseq$gene_type = paste0(tolower(substr(whole_nucseq$gene_names, 4, 4)), '_gene')
    whole_nucseq_by_gene = whole_nucseq[gene_type == GENE_NAME]
    sequence_composition = letterFrequency(DNAStringSet(whole_nucseq_by_gene$sequences), letters = c('T', 'G', 'C', 'A'))
    rownames(sequence_composition) = whole_nucseq_by_gene$gene_names
    frequencies = sequence_composition/rowSums(sequence_composition)
    return(frequencies[rownames(frequencies) == gene])
}

find_base_background_frequencies <- function(){
    #TODO not totally sure what to do here, but I think just find the
    #frequencies across ALL possible v-genes for example (but maybe I should be
    #doing this on a per-subject basis...)
    #OR compute an A/T and G/C background proportion...since pnucs deal with
    #either or...
    whole_nucseq = fread('_ignore/tcrb_processed_geneseq.tsv')
    whole_nucseq$gene_type = paste0(tolower(substr(whole_nucseq$gene_names, 4, 4)), '_gene')
    whole_nucseq_by_gene = whole_nucseq[gene_type == GENE_NAME]
    sequence_composition = letterFrequency(DNAStringSet(whole_nucseq_by_gene$sequences), letters = c('T', 'G', 'C', 'A'))
    background_frequency = colSums(sequence_composition)/sum(sequence_composition)
    return(background_frequency)
}

calculate_PPM_from_PFM <- function(PFM){
    PFM_filtered = PFM[c('T', 'G', 'C', 'A'),]
    PPM = t(t(PFM_filtered)/colSums(PFM_filtered))
    return(PPM)
}

calculate_PWM_from_PFM <- function(PFM, weighting = NULL){
    background_freqs = find_base_background_frequencies()
    background_matrix = matrix(background_freqs, length(background_freqs), (RIGHT_NUC_MOTIF_COUNT+LEFT_NUC_MOTIF_COUNT))

    PPM = calculate_PPM_from_PFM(PFM)
    rownames(background_matrix) = names(background_freqs)
    stopifnot(identical(rownames(background_matrix), rownames(PPM)))

    PWM = log2(PPM/background_matrix)
    return(PWM)
}

# old protocol: 
# PCM_all = compile_PFMs_old(group_localID_list = get_subject_partition_by_SNP(), group_name = get_group_name())
# PPM_all = calculate_PPM_from_PFM(PCM_all)
# PWM_all = calculate_PWM_from_PFM(PPM_all)

# old protocol one subject:
subject = "HIP00110"
PCM_all = compile_PFMs_old(group_localID_list = c(subject), group_name = 'test')
PPM_all = calculate_PPM_from_PFM(PCM_all)
PWM_all = calculate_PWM_from_PFM(PPM_all)

# new protocol one subject

PPM_all_new = compile_PPMs_new(group_localID_list = c(subject), group_name = 'test')

