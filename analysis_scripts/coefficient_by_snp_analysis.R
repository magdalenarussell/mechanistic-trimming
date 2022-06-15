source('config/config.R')

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
stopifnot(ANNOTATION_TYPE %in% c('igor', 'parsimony'))

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
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[8])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[9])

UPPER_TRIM_BOUND <<- as.numeric(args[10]) 

MODEL_TYPE <<- args[11]
stopifnot(MODEL_TYPE != 'null')

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[12])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/individual_comparison_functions.R')
source('analysis_scripts/snp_analysis_functions.R')

snps = c(20717849, 20717772)
snp_genotypes = compile_all_genotypes_snp_list(snps)  

#TODO remove this HIP conversion
gwas_id = fread(ID_MAPPING_FILE)
snp_genotypes = merge(snp_genotypes, gwas_id)[, -c('scanID', 'localID')]

# Compile data for all subjects
motif_data = aggregate_all_subject_data()
together = merge(motif_data, snp_genotypes, by = 'subject')
together = calculate_subject_gene_weight(together)

# Fit model
model = fit_model(together)

# get set of all variables
vars = names(coef(model))
vars_no_motif = vars[!(vars %like% 'motif')]
motif_terms = get_positions()
all_vars = c(vars_no_motif, motif_terms)
all_vars = str_replace_all(all_vars, '\\(', '\\\\(')
all_vars = str_replace_all(all_vars, '\\)', '\\\\)')

# vary each variable by snp and calculate model loss
for (snp in snps){
    new_colnames = str_replace(colnames(together), paste0(snp), 'temp_snp')
    snp_data = together
    colnames(snp_data) = new_colnames
    snp_data = snp_data[!is.na(temp_snp)]
    reset_colnames = str_replace(colnames(snp_data), 'temp_snp', paste0(snp))
    colnames(snp_data) = reset_colnames

    # redefine gene weighted with this new, subsetted group of subjects (only people with genotype data available)
    snp_data = calculate_subject_gene_weight(snp_data)

    # Fit model, write predicted distribution and pwm files
    normal_model = fit_model(snp_data)
    total_loglik = logLik(normal_model) * nrow(snp_data)

    snp_results = data.table()
    for (variable in all_vars){
        new_formula = include_snp_in_model_formula(get_model_formula, variable, snp)
        var_model = fit_model(snp_data, new_formula)        
        var_loglik = logLik(var_model) * nrow(snp_data)
        df_diff = length(coef(var_model)) - length(coef(normal_model))
        lrt_p = calculate_lrt_p(total_loglik, var_loglik, df_diff)
        temp_result = compile_snp_result(total_loglik, var_loglik, df_diff, lrt_p, variable, snp)
        snp_results = rbind(snp_results, temp_result)
    }
    write_snp_result_dt(snp_results, snp)
}
