LOCUS <<- 'delta'
source(paste0(MOD_PROJECT_PATH, '/scripts/locus_specific_functions/', LOCUS, '.R'))

extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]]
    len = length(file_name)
    file_name = file_name[len]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    return(file_root_name)
}

get_whole_nucseqs <- function(){
    whole_nucseq = fread(get(paste0('WHOLE_NUCSEQS_', 'delta')))[, -c('names')]
    return(whole_nucseq[, c('gene', 'sequence')])
}
