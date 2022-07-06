source('config/config.R')

library(foreach)
library(doParallel)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(devtools)
library(ggbiplot)
library(ggplot2)
library(cowplot)
library(mclogit)
library(matrixcalc)
library(RhpcBLASctl)
library(RColorBrewer)
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

MODEL_TYPE <<- args[11]

if (grepl('_side_terminal', MODEL_TYPE, fixed = TRUE) | grepl('two-side-base-count', MODEL_TYPE, fixed = TRUE) | grepl('left-base-count', MODEL_TYPE, fixed = TRUE)| grepl('two-side-dinuc-count', MODEL_TYPE, fixed = TRUE)){
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

coef_snps[!(base ==  ''), parameter := paste0(parameter, ', base ', base)]
coef_snps = coef_snps[!is.na(parameter)]

formatted_coef_snps = unique(coef_snps[, -c('base')]) %>% 
    pivot_wider(names_from = 'parameter', values_from = 'coefficient') %>%
    as.data.table()

file_path = get_individual_comparison_file_path()

for (snp in snps){
    snp_subset = formatted_coef_snps[!is.na(eval(as.name(snp)))]
    params = unique(coef_snps$parameter)
    pca = prcomp(snp_subset[,..params], center = TRUE,scale. = TRUE)
    
    color_palette = brewer.pal(3, 'Set2') 
    comparisons = list(c(1, 2), c(3, 4))
    for (comp in comparisons) {
        plot_name = file.path(file_path, paste0('snp', snp, '_pca_', comp,'.pdf'))
        plot_name2 = file.path(file_path, paste0('snp', snp, '_pca_', comp,'_no_points.pdf'))

        plot_pca = ggbiplot(pca, groups=as.factor(snp_subset[[paste(snp)]]), ellipse = TRUE, choices = comp, alpha = 0.7, labels.size = 4) +
            theme_cowplot(font_family = 'Arial') + 
            background_grid(major = 'xy') + 
            scale_color_manual(values = color_palette) +
            panel_border(color = 'gray60', size = 1.5) +
            scale_x_continuous(expand = expansion(mult = 0.4))+
            scale_y_continuous(expand = expansion(mult = 0.4))+
            theme(text = element_text(size = 20), axis.text = element_text(size = 15)) 
    
        plot_pca2 = ggbiplot(pca, groups=as.factor(snp_subset[[paste(snp)]]), ellipse = FALSE, choices = comp, alpha = 0, labels.size = 4) +
            theme_cowplot(font_family = 'Arial') + 
            background_grid(major = 'xy') + 
            scale_color_manual(values = color_palette) +
            panel_border(color = 'gray60', size = 1.5) +
            scale_x_continuous(expand = expansion(mult = 0.4))+
            scale_y_continuous(expand = expansion(mult = 0.4))+
            theme(text = element_text(size = 20), axis.text = element_text(size = 15)) 
 
        ggsave(plot_name, plot = plot_pca, width = 8, height = 8, units = 'in', dpi = 750, device = cairo_pdf)
        ggsave(plot_name2, plot = plot_pca2, width = 8, height = 8, units = 'in', dpi = 750, device = cairo_pdf)

    }
}
