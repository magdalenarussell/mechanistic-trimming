source('mechanistic-trimming/config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

ANNOTATION_TYPE <<- args[1]
stopifnot(ANNOTATION_TYPE %in% c('igor_alpha', 'igor_beta', 'adaptive_alpha', 'adaptive_beta'))

PARAM_GROUP <<- args[2]
param_types = list.files(path = paste0(MOD_PROJECT_PATH, '/scripts/param_groups/'))
param_types = str_sub(param_types, end = -3)
stopifnot(PARAM_GROUP %in% param_types)
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(args[3])

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[4])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[5])

MODEL_TYPE <<- args[6]

PROP <<- as.numeric(args[7])

source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))

# Read processed data
filename = processed_data_path()
motif_data = fread(filename)

registerDoParallel(cores=NCPU)
foreach(iter = seq(50)) %dopar% {
    temp = subsample(motif_data, PROP)
    name = subsampling_processed_data_path(PROP, iter)
    fwrite(temp, name, sep = '\t')
    print(paste0('finished with ', iter, ' of 50'))
}
stopImplicitCluster()

