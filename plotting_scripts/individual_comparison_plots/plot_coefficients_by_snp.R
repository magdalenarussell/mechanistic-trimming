source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(ggplot2)
library(cowplot)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)


args = commandArgs(trailingOnly=TRUE)

ANNOTATION_TYPE<<- args[1]
TRIM_TYPE <<- args[2]
trim_types = list.files(path = 'scripts/gene_specific_functions/')
trim_types = str_sub(trim_types, end = -3)
stopifnot(TRIM_TYPE %in% trim_types)

PRODUCTIVITY <<- args[3]

MOTIF_TYPE <<- args[4] 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[5])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- args[6]
GENE_WEIGHT_TYPE <<- args[7]

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[8])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[9])

UPPER_TRIM_BOUND <<- as.numeric(args[10]) 
LOWER_TRIM_BOUND <<- RIGHT_NUC_MOTIF_COUNT - 2 

MODEL_TYPE <<- args[11]

if (grepl('_side_terminal_melting', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[12])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/individual_comparison_functions.R')

# Read in model coefficient data 
data_dir = get_coefficient_output_file_path()
indiv_coefs = combine_all_individual_coefficients(data_dir)
# get snp genotypes for missense Artemis snp (20717849) and most significant artemis SNP for V- and J- gene trimming (20717772)
snps = c(20717849, 20717772)
snp_genotypes = compile_all_genotypes_snp_list(snps)  

#TODO remove this HIP conversion
gwas_id = fread(ID_MAPPING_FILE)
snp_genotypes = merge(snp_genotypes, gwas_id)[, -c('scanID', 'localID')]

coef_snps = merge(indiv_coefs, snp_genotypes, by.x = 'model_group', by.y = 'subject', all = TRUE)

formatted_coef_snps = coef_snps %>% 
    pivot_longer(starts_with("2"), names_to = 'snp', values_to = 'genotype') %>%
    as.data.table()
formatted_coef_snps[!(base ==  ''), parameter := paste0(parameter, ', base ', base)]

for (snp in snps) {
    for (parameter_group in c('motif', 'trim_length', 'terminal_melting')) {
        if (!(MODEL_TYPE %like% parameter_group)){
            if (parameter_group != 'trim_length' | !(MODEL_TYPE %like% 'distance')){
                next
            }
        }
        plot_coefficient_by_snp(formatted_coef_snps, snp, parameter_group)
    }
}

