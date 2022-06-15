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

PRODUCTIVITY <<- args[2]

MOTIF_TYPE <<- args[3] 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[4])

MODEL_GROUP <<- args[5]

GENE_WEIGHT_TYPE <<- args[6]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[7])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[8])

UPPER_TRIM_BOUND <<- as.numeric(args[9]) 

MODEL_TYPE <<- args[10]
stopifnot(MODEL_TYPE != 'null')

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[11])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

TRIM_TYPE <<- 'v_trim'
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('analysis_scripts/snp_analysis_functions.R')
source('analysis_scripts/v_j_comparison_functions.R')

# Compile data for all subjects
v_motif_data = aggregate_all_subject_data()

TRIM_TYPE <<- 'j_trim'
GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('analysis_scripts/snp_analysis_functions.R')

# Compile data for all subjects
j_motif_data = aggregate_all_subject_data()

v_motif_data$gene_type = 'v_gene'
j_motif_data$gene_type = 'j_gene'

together = rbind(v_motif_data, j_motif_data) 

# Fit model
model = fit_model(together)
original_loglik = logLik(model) * nrow(together)

# get set of all variables
vars = names(coef(model))
vars_no_motif = vars[!(vars %like% 'motif')]
motif_terms = get_positions()
all_vars = c(vars_no_motif, motif_terms)
all_vars = str_replace_all(all_vars, '\\(', '\\\\(')
all_vars = str_replace_all(all_vars, '\\)', '\\\\)')

all_var_results = data.table()
for (variable in all_vars){
    gene_variation_formula = include_gene_type_in_model_formula(get_model_formula, variable)
    gene_var_model = fit_model(together, gene_variation_formula)
    gene_var_loglik = logLik(gene_var_model) * nrow(together)
    
    df_diff = length(coef(gene_var_model)) - length(coef(model))
    lrt_p = calculate_lrt_p(original_loglik, gene_var_loglik, df_diff)
    
    temp_result = compile_gene_comparison_result(original_loglik, gene_var_loglik, df_diff, lrt_p, variable)

    all_var_results = rbind(all_var_results, temp_result)
}

write_gene_comparison_result_dt(all_var_results)
