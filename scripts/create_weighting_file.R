library(plyr)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)
library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)

extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][3]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_")[[1]][3]
    return(localID)
}

get_observed_counts <- function(file){
    data = fread(file)
    together = data.frame(matrix(ncol = 2, nrow = 0))
    colnames(together) = c('gene', extract_subject_ID(file))
    for (gene_type in c('v_gene', 'd_gene', 'j_gene')){
        data_N = data[,.N, by=gene_type]
        colnames(data_N) = c('gene', extract_subject_ID(file))
        together = rbind(together, data_N)
    }
    transformed = as.data.frame(t(together))
    colnames(transformed) = transformed[1,]
    transformed = mutate_all(transformed[-1,], function(x) as.numeric(as.character(x)))
    transformed$subject = extract_subject_ID(file)
    return(transformed)
}
    

create_all_PFMs <- function(directory, path){
    files = fs::dir_ls(path = directory)
    whole_nucseq = fread('_ignore/tcrb_processed_geneseq.tsv')
    colnames = whole_nucseq$gene_names[-c(1:24)]
    empirical_counts = data.frame(matrix(ncol = length(colnames), nrow = 0))
    colnames(empirical_counts) = colnames
    for (file in files){
        empirical_counts = rbind.fill(empirical_counts, get_observed_counts(file))
    }
    
    empirical_counts[is.na(empirical_counts)] = 0
    expected_empirical = summarise_all(empirical_counts, ~if(is.numeric(.)) sum(.) else "ALL") 
    empirical_counts = rbind(empirical_counts, expected_empirical)
    empirical_counts = empirical_counts %>% mutate(total_tcr = rowSums(across(where(is.numeric))))
    empirical_frequencies = empirical_counts %>% mutate_at(vars(-c('subject', 'total_tcr')), ~(./total_tcr))
    write.table(empirical_counts, paste0(path, 'empirical_gene_counts.tsv'), sep = '\t')
    write.table(empirical_frequencies, paste0(path, 'empirical_gene_frequencies.tsv'), sep = '\t')
}

create_all_PFMs(directory = "_ignore/emerson_stats", path = '/fh/fast/matsen_e/shared/tcr-gwas/exomotif/')
