source('mechanistic-trimming/config/config.R')

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

ANNOTATION_TYPE <<- 'igor_sim_alpha' 

PARAM_GROUP <<- 'nonproductive_v-j_trim'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

MODEL_TYPE <<- 'motif_two-side-base-count-beyond_mh'

L2 <<- 'True'

source(paste0(MOD_PROJECT_PATH, '/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/model_fitting_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

filename = processed_data_path()
motif_data = fread(filename)

plot_mh_by_trim_length <- function(motif_data, trim_type, standardize = FALSE, center = FALSE){
    cols = colnames(motif_data)[colnames(motif_data) %like% 'mh_prop']
    cols = c(trim_type, cols)
    subset = motif_data[, ..cols]

    long_subset = subset %>% pivot_longer(starts_with('mh_prop'), names_to = "mh_position", values_to = 'mh') %>% as.data.table()
    long_subset[, overlap := substring(mh_position, nchar(mh_position), nchar(mh_position))]
    long_subset[mh_position %like% 'up', position := 'up']
    long_subset[mh_position %like% 'mid', position := 'mid']
    long_subset[mh_position %like% 'down', position := 'down']
    long_subset$position = factor(long_subset$position, levels = c('up', 'mid', 'down'))

    mh_var = 'mh'
    if (isTRUE(standardize)){
        long_subset[, mean := mean(mh)]
        long_subset[, sd := sd(mh)]
        long_subset[, std_mh := (mh - mean)/sd]
        mh_var = 'std_mh'
    }

    if (isTRUE(center)){
        long_subset[, mean := mean(mh)]
        long_subset[, centered_mh := (mh - mean)]
        mh_var = 'centered_mh'
    }

    long_subset_sample = long_subset[sample(.N, 100000)]

    plot = ggplot(long_subset_sample) +
        geom_point(aes(x = get(trim_type), y = get(mh_var), color = overlap), size = 4, alpha = 0.5) +
        facet_wrap(~position)+
        theme_cowplot(font_family = 'Arial') + 
        theme(text = element_text(size = 35), axis.text.x=element_text(size = 25), axis.text.y = element_text(size = 25), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) +
        xlab(trim_type) +
        ylab(mh_var)
    return(plot)
}

# for V-trim
plot = plot_mh_by_trim_length(motif_data, 'v_trim')
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_by_trim_length/v_trim.pdf')

ggsave(file_name, plot = plot, width = 15, height = 6, units = 'in', dpi = 750, device = cairo_pdf)

# for J-trim
plot = plot_mh_by_trim_length(motif_data, 'j_trim')
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_by_trim_length/j_trim.pdf')

ggsave(file_name, plot = plot, width = 15, height = 6, units = 'in', dpi = 750, device = cairo_pdf)

# for V-trim
plot = plot_mh_by_trim_length(motif_data, 'v_trim', standardize = TRUE)
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_by_trim_length/v_trim_std.pdf')

ggsave(file_name, plot = plot, width = 15, height = 6, units = 'in', dpi = 750, device = cairo_pdf)

# for J-trim
plot = plot_mh_by_trim_length(motif_data, 'j_trim', standardize = TRUE)
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_by_trim_length/j_trim_std.pdf')

ggsave(file_name, plot = plot, width = 15, height = 6, units = 'in', dpi = 750, device = cairo_pdf)

# for V-trim
plot = plot_mh_by_trim_length(motif_data, 'v_trim', center = TRUE)
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_by_trim_length/v_trim_center.pdf')

ggsave(file_name, plot = plot, width = 15, height = 6, units = 'in', dpi = 750, device = cairo_pdf)

# for J-trim
plot = plot_mh_by_trim_length(motif_data, 'j_trim', center = TRUE)
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_by_trim_length/j_trim_center.pdf')

ggsave(file_name, plot = plot, width = 15, height = 6, units = 'in', dpi = 750, device = cairo_pdf)
