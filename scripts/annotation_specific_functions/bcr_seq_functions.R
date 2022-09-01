extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][3]
    localID = str_split(file_name, "_db.tab")[[1]][1]
    return(localID)
}

get_whole_nucseqs <- function(){
    whole_nucseq = fread(get('WHOLE_NUCSEQS_igh'))[, -c('name')]
    return(whole_nucseq[, c('gene', 'sequence')])
}
