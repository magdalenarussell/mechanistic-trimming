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

PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

MODEL_TYPE <<- 'motif_two-side-base-count-beyond_ligation-mh'

L2 <<- 'False'

LIGATION_PARAM <<- 0
TRIMMING_PROB_MODEL <<- 'motif_two-side-base-count-beyond'
ANNOTATION_SUBTYPE <<- 'always-ligate'

ANNOTATION_TYPE <<- 'igor_mh_sim_alpha'
old_annotation = ANNOTATION_TYPE

source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

ANNOTATION_TYPE <<- paste0(old_annotation, '_from_', TRIMMING_PROB_MODEL, '_MHprob', LIGATION_PARAM , '_', ANNOTATION_SUBTYPE)
# Read in model coefficient data 
coef_path = get_model_coef_file_path(L2)
coefs = fread(coef_path)

# get single V params
PARAM_GROUP <<- 'nonproductive_v_trim'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

MODEL_TYPE <<- 'motif_two-side-base-count-beyond'

L2 <<- 'False'

ANNOTATION_TYPE <<- 'igor_sim_alpha'

source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

# Read in model coefficient data 
coef_path = get_model_coef_file_path(L2)
Vcoefs = fread(coef_path)

# get single J params
PARAM_GROUP <<- 'nonproductive_j_trim'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

source(paste0(MOD_PROJECT_PATH, '/config/file_paths.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

# Read in model coefficient data 
coef_path = get_model_coef_file_path(L2)
Jcoefs = fread(coef_path)

#widen data
coefs$model = "Joint"
Vcoefs$model = "Single"
Jcoefs$model = "Single"

coefs = coefs %>% pivot_wider(names_from = model, values_from = value) %>% as.data.table()
Vcoefs = Vcoefs %>% pivot_wider(names_from = model, values_from = value) %>% as.data.table()
Jcoefs = Jcoefs %>% pivot_wider(names_from = model, values_from = value) %>% as.data.table()

single = rbind(Vcoefs, Jcoefs)
tog = merge(coefs, single, by = c('coefficient', 'base', 'position', 'side', 'trim_type'))

fwrite(coefs, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/compare_single_gene_coef/coefs_MHprob0_single_gene_comparison.tsv'), sep = '\t')

tog$type = paste0(tog$side, ', ', tog$coefficient, ', ', tog$base)

scatter = ggplot(tog)+
    geom_abline(intercept = 0, slope = 1, size = 3, color = 'gray60', linetype = 'dashed') +
    geom_point(aes(x = Joint/log(10), y = Single/log(10), color = type, shape = trim_type), size = 10) +
    xlab('log10(trimming probability from joint/MH model)')+
    ylab('log10(trimming probability from single gene models)')+
    theme_cowplot(font_family = 'Arial') + 
    labs(color = 'Coefficient type', shape = 'Trimming type') +
    background_grid(major = 'xy') + 
    # scale_color_brewer(palette = 'Set2') + 
    panel_border(color = 'gray60', size = 1.1) +
    theme(text = element_text(size = 18), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 16), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), panel.spacing = unit(3, "lines"), legend.key.height = unit(2, 'cm'), legend.key.width = unit(1.5, 'cm'))  

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/mh_simulator_experiments/compare_single_gene_coef/scatter_MHprob0_single_gene_comparison.pdf')

ggsave(file_name, plot = scatter, width = 18, height = 14, units = 'in', dpi = 750, device = cairo_pdf)

 
