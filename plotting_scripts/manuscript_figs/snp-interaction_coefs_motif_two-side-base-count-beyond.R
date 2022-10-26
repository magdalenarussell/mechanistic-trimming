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
library(ggseqlogo)
omp_set_num_threads(1)
blas_set_num_threads(1)

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
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

UPPER_TRIM_BOUND <<- as.numeric(14) 

MODEL_TYPE <<- 'motif_two-side-base-count-beyond_snp-interaction-20717772'

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- 10

source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')
source('plotting_scripts/individual_comparison_functions.R')

filename =  get_model_bootstrap_file_name() 
pwm = fread(filename)
pwm = pwm[parameter %like% 'snp']
pwm$parameter = str_replace(pwm$parameter, ':snp', '')
pwm$parameter = str_replace(pwm$parameter, '_prop', '')
pwm$coefficient = -1*(pwm$coefficient)

# plot_model_coefficient_heatmap(pwm, with_values = TRUE, limits = c(-0.4, 0.4))
heatmap = plot_model_coefficient_heatmap_single_group(pwm, with_values = FALSE, write_plot = FALSE, limits = c(-0.011, 0.011))
heatmap = heatmap + 
    theme(legend.position = 'none', text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20)) 

heatmap2 = plot_base_count_coefficient_heatmap_single_group(pwm, with_values = FALSE, write_plot = FALSE, limits = c(-0.011, 0.011))
heatmap2 = heatmap2 + 
    theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 20),legend.position = 'bottom', legend.direction = 'horizontal', legend.justification="center") +
    guides(fill = guide_colourbar(barwidth = 27.5, barheight = 2)) +
    labs(fill = 'log10(probability of deletion)\t')

legend = get_legend(heatmap2) 
heatmap2 = heatmap2 + theme(legend.position = 'none')

path = get_manuscript_path()
file_name = paste0(path, '/snp_interaction_heatmaps.pdf')
ggsave(file_name, plot = first_grid, width = 20, height = 5.5, units = 'in', dpi = 750, device = cairo_pdf)

# make logo
source('plotting_scripts/logo_functions.R')

ppm = convert_pwm_to_ppm(pwm)
ppm_matrix = convert_ppm_to_matrix(ppm)
position_values = map_positions_to_values(unique(ppm$parameter))

logo = ggplot() + 
    # geom_logo(ppm_matrix, col_scheme='base_pairing') + 
    geom_logo(ppm_matrix) + 
    theme_logo() +    
    geom_vline(xintercept = 1.5, size = 4, color = 'black')+
    annotate("text", x = 0.5, y = 0, label = "5\'", size = 8) +  
    annotate("text", x = 3.6, y = 0, label = "3\'", size = 8) +  
    theme_cowplot(font_family = 'Arial') + 
    xlab('Position') +
    ylab ('Bits') +
    theme(text = element_text(size = 30), legend.position = 'none', axis.text = element_text(size = 20), axis.line = element_line(color = 'gray60', size = 1.5), axis.ticks = element_line(color = 'gray60', size = 0.75)) +
    coord_cartesian(ylim = c(0, 1.6e-05), clip = "off")

    # panel_border(color = 'gray60', size = 1.5) 

logo$scales$scales[[1]] <- scale_x_continuous(breaks= c(1, 2, 3),labels=c("-1", "1", "2"))

all = align_plots(heatmap, heatmap2, logo, legend, align = 'vh', axis = 'lbr')

first_grid = plot_grid(all[[1]], all[[2]], nrow = 1, labels = c("B", "C"), label_size = 35, rel_widths = c(1, 1), align = 'h')
first_grid = plot_grid(first_grid, NULL, legend, nrow = 3, rel_heights = c(1, 0.05,0.2)) 
all_logo = plot_grid(all[[3]], NULL, nrow = 2, rel_heights = c(1, 0.25))
tog = plot_grid(NULL, all_logo, NULL, first_grid, labels = c('A', '', '', ''), rel_widths = c(0.015,0.3, 0.005,0.6), label_size = 35, nrow = 1)

path = get_manuscript_path()
file_name = paste0(path, '/snp_interaction_heatmaps.pdf')
ggsave(file_name, plot = tog, width = 20, height = 6, units = 'in', dpi = 750, device = cairo_pdf)


