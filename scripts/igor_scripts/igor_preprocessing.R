library(data.table)
library(tidyverse)

source('config/config.R')
source('scripts/igor_scripts/igor_processing_functions.R')

args = commandArgs(trailingOnly=TRUE)

file = args[1]
output_directory = args[2]

subject_ID = extract_subject_ID(file)
cdr3s = get_raw_cdr3_seqs(file)

output_location = file.path(output_directory, subject_ID)
dir.create(output_location, recursive = TRUE, showWarnings = FALSE)

fwrite(as.data.table(cdr3s), file.path(output_location, 'cdr3_seqs.txt'), sep = '\t', col.names = FALSE)

cat(output_location)
