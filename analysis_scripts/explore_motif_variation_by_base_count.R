source('config/config.R')

library(ggplot2)
library(cowplot)
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
stopifnot(TRIM_TYPE == 'v_trim')

PRODUCTIVITY <<- args[3]

MOTIF_TYPE <<- args[4] 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(args[5])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- args[6]

GENE_WEIGHT_TYPE <<- args[7]
stopifnot(GENE_WEIGHT_TYPE %in% c('p_gene_given_subject', 'p_gene_marginal', 'raw_count', 'uniform'))

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[8])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[9])

UPPER_TRIM_BOUND <<- as.numeric(args[10]) 

MODEL_TYPE <<- args[11]
stopifnot(MODEL_TYPE %like% 'motif')
stopifnot(MODEL_TYPE %like% 'base-count')

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- as.numeric(args[12])
} else {
    LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA
}

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')
source('analysis_scripts/pwm_profile_functions.R')

# Compile data for all subjects
motif_data = aggregate_all_subject_data()

# find genes with the same base count profiles
positions = get_positions()
left_vars = get_all_base_variables('left')
right_vars = get_all_base_variables('right')
cols = c('gene', 'p_gene', 'trim_length', 'motif', positions, left_vars, right_vars)
condensed = unique(motif_data[,..cols]) 

# Read in model coefficient data 
pwm = get_model_coefficient_data() 

# calculate pwm score by gene and trim length
condensed[, pwm_score := unlist(lapply(motif, function(x) as.numeric(get_pwm_score(pwm, x, positions))))]

small_diffs = data.table()
large_diffs = data.table()
gene1_cols = c('gene', 'p_gene', 'trim_length', 'motif', positions, left_vars, right_vars, 'pwm_score') 
gene2_cols = c('gene', 'trim_length', 'motif', positions, 'pwm_score') 

for (gene1 in unique(condensed$gene)){
    for (gene2 in unique(condensed$gene)){
        temp1 = condensed[gene == gene1][order(trim_length)]
        temp2 = condensed[gene == gene2][order(trim_length)]
        # make sure that the base count values are all equal
        if (all(temp1$right_base_count_AT == temp2$right_base_count_AT) & all(temp1$left_base_count_AT == temp2$left_base_count_AT) & !(all(temp1$pwm_score == temp2$pwm_score))) {
            if ((gene1 != gene2)){
                temp_1 = condensed[gene == gene1][, ..gene1_cols]
                temp_2 = condensed[gene == gene2][, ..gene2_cols]
                setnames(temp_2, 'gene', 'similar_gene')
                setnames(temp_2, 'motif', 'similar_gene_motif')
                setnames(temp_2, 'pwm_score', 'similar_gene_pwm_score')
                setnames(temp_2, positions, paste0('similar_gene_',positions))
                temp_2$gene = gene1
                temp = merge(temp_1, temp_2, by = c('gene', 'trim_length')) 
                # make sure that not all of the pwm score values (for trims greater than 3) are equal
                if (!(all(temp1[trim_length > 3]$pwm_score == temp2[trim_length > 3]$pwm_score))){
                    if (nrow(large_diffs) != 0){
                        if (nrow(large_diffs[gene == gene2 & similar_gene == gene1]) > 0) {
                            temp = data.table()
                        }
                    }
                    large_diffs = rbind(large_diffs, temp)                        
                } else {
                    if (nrow(small_diffs) != 0){
                        if (nrow(small_diffs[gene == gene2 & similar_gene == gene1]) > 0) {
                            temp = data.table()
                        }
                    }

                    small_diffs = rbind(small_diffs, temp)
                }
            }
        }
    }
}

path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, MODEL_GROUP, GENE_WEIGHT_TYPE, MODEL_TYPE , paste0(TRIM_TYPE, '_', MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
dir.create(path, recursive = TRUE)
 
############# explore how variations in motif when base count parameters remain constant
diffs = small_diffs[motif != similar_gene_motif]
diffs[, pwm_score_shift := similar_gene_pwm_score- pwm_score]



############# plot right GC count and left GC count by trim length 

condensed[, right_base_count_GC_by_trim := right_base_count_GC/trim_length]
plot = ggplot(condensed[order(trim_length)]) +
    geom_point(aes(x = trim_length, y = right_base_count_GC_by_trim), size = 3, alpha = 0.5) +
    geom_smooth(method = 'loess', aes(x = trim_length, y = right_base_count_GC_by_trim, group = gene), size = 3, alpha = 0.5, color = 'black', se = FALSE) +
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

file_name = file.path(path, 'compare_right_GC_count_by_length.pdf') 
ggsave(file_name, plot = plot, width = 22, height = 6, units = 'in', dpi = 750, device = cairo_pdf)

condensed[, left_base_count_GC_by_trim := left_base_count_GC/LEFT_SIDE_TERMINAL_MELT_LENGTH]
plot = ggplot(condensed[order(trim_length)]) +
    geom_point(aes(x = trim_length, y = left_base_count_GC_by_trim), size = 3, alpha = 0.5) +
    geom_smooth(method = 'loess', aes(x = trim_length, y = left_base_count_GC_by_trim, group = gene), size = 3, alpha = 0.5, color = 'black', se = FALSE) +
    theme_cowplot(font_family = 'Arial') + 
    theme(text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) 

file_name = file.path(path, 'compare_left_GC_count_by_length.pdf') 
ggsave(file_name, plot = plot, width = 22, height = 6, units = 'in', dpi = 750, device = cairo_pdf)


