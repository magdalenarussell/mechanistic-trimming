library(data.table)
library(doParallel)
library(foreach)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
gamma_dir = args[1]
output_dir = args[1]
dir.create(output_dir)

process_file <- function(file){
    file_data = fread(file)
    cols = c('v_resolved', 'd_resolved', 'j_resolved', 'v_deletions', 'd5_deletions', 'd3_deletions', 'j_deletions', 'n1_insertions', 'n2_insertions', 'frame_type')
    file_data = file_data[,..cols]
    new_names = c('v_gene', 'd_gene', 'j_gene', 'v_trim', 'd0_trim', 'd1_trim', 'j_trim', 'vd_insert', 'dj_insert', 'frame_type')
    colnames(file_data) = new_names
    file_data[frame_type == 'In', productive := TRUE]
    file_data[frame_type != 'In', productive := FALSE]
    file_data = file_data[v_gene != '']
    file_data[, v_gene := str_replace(v_gene, 'TCRG', 'TRG')]
    file_data[, j_gene := str_replace(j_gene, 'TCRG', 'TRG')]
    file_data[, d_gene := str_replace(d_gene, 'TCRG', 'TRG')]
    # if no allele was inferred, then assume allele 1
    file_data[!(v_gene %like% '\\*'), v_gene := paste0(v_gene, '*01')]
    file_data[!(j_gene %like% '\\*'), j_gene := paste0(j_gene, '*01')]
    return(file_data[, -c('frame_type')])
}

files = fs::dir_ls(path = gamma_dir)
registerDoParallel(cores=NCPU)
foreach(file = files) %dopar% {
    print(paste(file))
    file_name = str_split(file, '/')[[1]][8]
    processed = process_file(file) 
    fwrite(processed, file.path(output_dir, file_name), sep = '\t')
}
stopImplicitCluster()
 
