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
stopifnot(TRIM_TYPE == 'v_trim')

MOTIF_TYPE <<- args[3] 
stopifnot(MOTIF_TYPE %in% c('bounded', 'unbounded', 'unbounded_no_pnuc'))

NCPU <<- as.numeric(args[4])

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

MODEL_GROUP <<- args[5]
GENE_WEIGHT_TYPE <<- args[6]

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(args[7])
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(args[8])

UPPER_TRIM_BOUND <<- as.numeric(args[9]) 
LOWER_TRIM_BOUND <<- RIGHT_NUC_MOTIF_COUNT - 2 

MODEL_TYPE <<- 'motif' 

LEFT_SIDE_TERMINAL_MELT_LENGTH <<- NA


source('scripts/data_compilation_functions.R')
source('scripts/model_fitting_functions.R')
source('plotting_scripts/plotting_functions.R')

# Compile data for all subjects
motif_data = aggregate_all_subject_data()

positions = get_positions()
composition = data.table(base = c('A', 'T', 'C', 'G'))

for (pos in positions){
    temp = motif_data[, sum(count), by = get(pos)]
    colnames(temp) = c('base', pos)
    composition = merge(composition, temp)
}

data = composition %>%
        pivot_longer(!c(base), names_to = 'position', values_to = 'count') %>%
        as.data.table()

data[, freq := count/sum(count), by = position]
data[, overall_base_count := sum(count), by = base]
data[, overall_base_freq := overall_base_count/sum(overall_base_count), by = position]
data[, pwm := log2(freq/overall_base_freq)]
position_values = map_positions_to_values(unique(data$position))
together = merge(data, position_values, by.x = 'position', by.y = 'positions')

# order variables
together$base = factor(together$base, levels = c('T', 'G', 'C', 'A'))
together$values = factor(together$values)

motif_length = RIGHT_NUC_MOTIF_COUNT + LEFT_NUC_MOTIF_COUNT

plot = ggplot(together, aes(x=values, y=base, fill=pwm)) +
       geom_tile() +
       theme_cowplot(font_family = 'Arial') + 
       xlab('Position') +
       ylab ('Base') +
       theme(text = element_text(size = 20), axis.line = element_blank(), axis.ticks = element_blank()) +
       geom_vline(xintercept = LEFT_NUC_MOTIF_COUNT + 0.5, size = 3, color = 'black') +
       guides(fill = guide_colourbar(barheight = 14)) +
       scale_fill_distiller(palette = 'PuOr', name = 'log2(position base frequency/\noverall base frequency)', limits = c(-1.6, 1.6)) +
       annotate("text", x = 0.35, y = 0.25, label = "5\'", size = 8) +  
       annotate("text", x = motif_length + 0.65, y = 0.25, label = "3\'", size = 8) +  
       coord_cartesian(ylim = c(1, 4), clip = "off") + 
       geom_text(data = together, aes(x = values, y = base, label = round(pwm, 3)))

path = file.path(PROJECT_PATH, 'plots', ANNOTATION_TYPE, MODEL_GROUP, 'base_frequency', paste0(MOTIF_TYPE, '_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND))
dir.create(path, recursive = TRUE)

ggsave(paste0(path, '/freq_heatmap.pdf'), plot = plot, width = 9, height = 4.2, units = 'in', dpi = 750, device = cairo_pdf)


