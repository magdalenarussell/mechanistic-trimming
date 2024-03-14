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

PARAM_GROUP <<- 'nonproductive_v-j_trim_ligation-mh'
source(paste0(MOD_PROJECT_PATH, '/scripts/param_groups/', PARAM_GROUP, '.R'))

NCPU <<- as.numeric(2)

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- as.numeric(1)
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- as.numeric(2)

MODEL_TYPE <<- 'motif_two-side-base-count-beyond_ligation-mh'

L2 <<- 'False'

source(paste0(MOD_PROJECT_PATH,'/analysis_scripts/ligation-mh_signal_simulator_scripts/ligation-mh_simulator_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/scripts/data_compilation_functions.R'))
source(paste0(MOD_PROJECT_PATH,'/plotting_scripts/plotting_functions.R'))

# get all possible configurations
configs = read_frames_data()

# get possible configs after MH adjustment
zero_configs = configs[ligation_mh == 0]

# adjust configs for MH
converted_configs = get_possible_mh(zero_configs, keep_gene_seqs = FALSE)
converted_configs = count_mh_bordering_trim(converted_configs)
converted_configs = reassign_trimming_sites_with_mh(converted_configs)

cols = c('v_gene', 'j_gene', 'v_trim', 'j_trim', 'adjusted_v_trim', 'adjusted_j_trim', 'ligation_mh', 'frame_type', 'frame_stop')
converted_subset = converted_configs[, ..cols]


# compile dataframe
## get zero/nonzero MH total possible configs
configs_counts = configs[, .N, by = ligation_mh]
configs_counts[ligation_mh == 0, zero_nonzero := 'zero MH']
configs_counts[ligation_mh > 0, zero_nonzero := 'nonzero MH']
df = configs_counts[, sum(N), by = zero_nonzero]

df[, type := "Total possible configurations"]

## get zero/nonzero MH total possible annotations
annot_counts = zero_configs[, .N, by = ligation_mh]
annot_counts[ligation_mh == 0, zero_nonzero := 'zero MH']
annot_counts[ligation_mh > 0, zero_nonzero := 'nonzero MH']
annot_counts = annot_counts[, sum(N), by = zero_nonzero]
if (nrow(annot_counts[zero_nonzero == 'nonzero MH']) == 0){
    temp = data.table(zero_nonzero = 'nonzero MH', V1 = 0)
    annot_counts = rbind(annot_counts, temp)
}

annot_counts[, type := 'Total possible annotations\nbefore MH adjustment']
df = rbind(df, annot_counts)

## get zero/nonzero MH total possible annotations after MH adjustment
converted_condensed = converted_subset[, .N, by = .(v_gene, j_gene, adjusted_v_trim, adjusted_j_trim, ligation_mh)]

converted_counts = converted_condensed[, .N, by = ligation_mh]
converted_counts[ligation_mh == 0, zero_nonzero := 'zero MH']
converted_counts[ligation_mh > 0, zero_nonzero := 'nonzero MH']
converted_counts = converted_counts[, sum(N), by = zero_nonzero]

converted_counts[, type := 'Total possible annotations\nafter MH adjustment']
df = rbind(df, converted_counts)

setnames(df, 'V1', 'count')

df[, type := factor(type, levels = c('Total possible configurations', 'Total possible annotations\nbefore MH adjustment', 'Total possible annotations\nafter MH adjustment'))]

plot = ggplot(df) +
    geom_bar(aes(x = type, y = count, fill = zero_nonzero), stat = 'identity')+
    ylab("Trimming/Ligation configuration count") +
    xlab("") +
    theme_cowplot(font_family = 'Arial') + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_brewer(palette = 'Set2') +
    theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 24))+
    labs(fill = 'Ligation configuration')

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/config_conversion.pdf')

ggsave(file_name, plot = plot, width = 20, height = 10, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE) 

fwrite(df, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/config_counts.tsv'), sep = '\t') 


## same plot, but for only nonprods
# compile dataframe
## get zero/nonzero MH total possible configs
configs_counts = configs[frame_type == 'Out' | frame_stop == TRUE, .N, by = ligation_mh]
configs_counts[ligation_mh == 0, zero_nonzero := 'zero MH']
configs_counts[ligation_mh > 0, zero_nonzero := 'nonzero MH']
df = configs_counts[, sum(N), by = zero_nonzero]

df[, type := "Total possible configurations"]

## get zero/nonzero MH total possible annotations
annot_counts = zero_configs[frame_type == 'Out' | frame_stop == TRUE, .N, by = ligation_mh]
annot_counts[ligation_mh == 0, zero_nonzero := 'zero MH']
annot_counts[ligation_mh > 0, zero_nonzero := 'nonzero MH']
annot_counts = annot_counts[, sum(N), by = zero_nonzero]
if (nrow(annot_counts[zero_nonzero == 'nonzero MH']) == 0){
    temp = data.table(zero_nonzero = 'nonzero MH', V1 = 0)
    annot_counts = rbind(annot_counts, temp)
}

annot_counts[, type := 'Total possible annotations\nbefore MH adjustment']
df = rbind(df, annot_counts)

## get zero/nonzero MH total possible annotations after MH adjustment
converted_condensed = converted_subset[frame_type == 'Out' | frame_stop == TRUE, .N, by = .(v_gene, j_gene, adjusted_v_trim, adjusted_j_trim, ligation_mh)]

converted_counts = converted_condensed[, .N, by = ligation_mh]
converted_counts[ligation_mh == 0, zero_nonzero := 'zero MH']
converted_counts[ligation_mh > 0, zero_nonzero := 'nonzero MH']
converted_counts = converted_counts[, sum(N), by = zero_nonzero]

converted_counts[, type := 'Total possible annotations\nafter MH adjustment']
df = rbind(df, converted_counts)

setnames(df, 'V1', 'count')

df[, type := factor(type, levels = c('Total possible configurations', 'Total possible annotations\nbefore MH adjustment', 'Total possible annotations\nafter MH adjustment'))]

plot = ggplot(df) +
    geom_bar(aes(x = type, y = count, fill = zero_nonzero), stat = 'identity')+
    ylab("Trimming/Ligation configuration count") +
    xlab("") +
    theme_cowplot(font_family = 'Arial') + 
    background_grid(major = 'xy') + 
    panel_border(color = 'gray60', size = 1.5) +
    scale_fill_brewer(palette = 'Set2') +
    theme(text = element_text(size = 30), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 24))+
    labs(fill = 'Ligation configuration')

# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/config_conversion_nonprod.pdf')

ggsave(file_name, plot = plot, width = 20, height = 10, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE) 

fwrite(df, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/config_counts_nonprod.tsv'), sep = '\t') 


# plot converted sites
converted_sites = converted_subset[ligation_mh > 0]
gene_choice_probs = get_igor_gene_usage_params()
top_V = gene_choice_probs$v_choice[order(-v_gene_prob)][1:6]
top_J = gene_choice_probs$j_choice[order(-j_gene_prob)][1:6]

converted_df = converted_sites[v_gene %in% top_V$v_gene & j_gene %in% top_J$j_gene, .N, by = .(v_gene, j_gene, v_trim, j_trim)]
converted_df[, ligation_mh := 0]

all_options_plot = ggplot(converted_df) +
                   facet_grid(cols = vars(j_gene), rows = vars(v_gene))+
                   geom_tile(aes(x = v_trim, y = j_trim), fill = 'black') +
                   geom_point(aes(x = v_trim, y = j_trim, color = ligation_mh)) +
                   ylab("J-gene trimming length (before MH adjustment)") +
                   xlab("V-gene trimming length (before MH adjustment)") +
                   theme_cowplot(font_family = 'Arial') + 
                   background_grid(major = 'xy') + 
                   panel_border(color = 'gray60', size = 1.5) +
                   theme(text = element_text(size = 40), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_text(size = 30), legend.key.size = unit(3, "cm"))+
                   labs(fill = 'Config with MH\nadjustment option', color = 'Zero MH\nligation option')


# save plot
file_name = paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/converted_sites_top_genes.pdf')

ggsave(file_name, plot = all_options_plot, width = 35, height = 30, units = 'in', dpi = 750, device = cairo_pdf, limitsize=FALSE) 

# save data               
fwrite(converted_df, paste0(MOD_PROJECT_PATH, '/plotting_scripts/manuscript_figs/ligation_mh_choices_germline_tra/converted_sites_top_genes.tsv'), sep = '\t') 


