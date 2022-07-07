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

ANNOTATION_TYPE<<- 'igor' 
TRIM_TYPE <<- 'v_trim' 
trim_types = list.files(path = 'scripts/gene_specific_functions/')
trim_types = str_sub(trim_types, end = -3)
stopifnot(TRIM_TYPE %in% trim_types)

PRODUCTIVITY <<- 'nonproductive' 

MOTIF_TYPE <<- 'unbounded' 
motif_types = list.files(path = 'scripts/motif_class_functions/')
motif_types = str_sub(motif_types, end = -3)
stopifnot(MOTIF_TYPE %in% motif_types)

NCPU <<- as.numeric(2)

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')

MODEL_GROUP <<- 'all_subjects' 
GENE_WEIGHT_TYPE <<- 'p_gene_marginal' 

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(2)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(4)

UPPER_TRIM_BOUND <<- as.numeric(14) 

MODEL_TYPE <<- 'motif' 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')

# get motif data
motif_data = aggregate_all_subject_data()
ordered_genes = motif_data[, .N, by = .(gene, p_gene)][order(-p_gene)]
top_genes = ordered_genes[1:6]$gene

# Read in dist data
predicted_trims = get_predicted_distribution_data() 

for (gene in unique(top_genes)){
    index = which(top_genes == gene)
    temp_plot = plot_predicted_trimming_dists(predicted_trims, gene, ylim = 0.75)
    temp_plot = temp_plot + 
        ggtitle(gene) +
        xlab('') +
        ylab('') +
        ylim(0, 0.8)
    assign(paste0('gene', index), temp_plot)
}

# Read in model coefficient data 
pwm = get_model_coefficient_data() 

# plot_model_coefficient_heatmap(pwm, with_values = TRUE, limits = c(-0.4, 0.4))
heatmap = plot_model_coefficient_heatmap_single_group(pwm, with_values = FALSE, write_plot = FALSE)
heatmap = heatmap + 
    theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20)) 


all = align_plots(gene1, gene2, gene3, gene4, gene5, gene6, heatmap, align = 'v', axis = 'l')
first_grid = plot_grid(all[[1]], all[[2]], all[[3]], all[[4]], all[[5]], all[[6]], nrow = 2) +
    draw_label("Number of trimmed nucleotides", x=0.5, y=  0, vjust=0 , angle= 0, size = 30, fontfamily = 'Arial') +
    draw_label("Probability", x=  0, y=0.5, vjust= 0.5, angle=90, size = 30, fontfamily = 'Arial')

temp_grid = plot_grid(NULL, first_grid, ncol = 2, rel_widths = c(0.017, 1))
second_grid = plot_grid(temp_grid, NULL, all[[7]], nrow = 3, labels = c("A", "", "B"), label_size = 35, rel_heights = c(1, 0.05, 0.7), align = 'v')

path = get_manuscript_path()
file_name = paste0(path, '/murugan_compare.pdf')
ggsave(file_name, plot = second_grid, width = 20, height = 14, units = 'in', dpi = 750, device = cairo_pdf)


