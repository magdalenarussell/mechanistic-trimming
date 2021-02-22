library(ggplot2)
library(tidyverse)
library(data.table)
setDTthreads(1)
library(Biostrings)
library(gdsfmt)
library(SNPRelate)
library(GWASTools)
library(RhpcBLASctl)
omp_set_num_threads(1)
blas_set_num_threads(1)

args = commandArgs(trailingOnly=TRUE)

TRIM_TYPE <<- args[1]
stopifnot(TRIM_TYPE == 'v_trim')

MOTIF_TYPE <<- args[2] 
stopifnot(MOTIF_TYPE %in% c('pnuc_motif', 'no_pnuc_motif'))

PFM_TYPE <<- args[3]
stopifnot(PFM_TYPE %in% c('unbounded'))

NCPU <<- args[4]

GENE_NAME <<- paste0(substring(TRIM_TYPE, 1, 1), '_gene')
stopifnot(GENE_NAME == 'v_gene')

# 5' motif nucleotide count
LEFT_NUC_MOTIF_COUNT <<- 4
# 3' motif nucleotide count
RIGHT_NUC_MOTIF_COUNT <<- 4

source('scripts/pfm_functions.R')
source('scripts/plot_functions.R')


plot_single_PWM_heatmap(group_name = get_group_name(), weighting = NULL)

significant_artemis_v_trim = compile_snps_from_GWAS(phenotype_list = c('v_trim'), gene = 'artemis')
ordered_top10_snps = significant_artemis_v_trim[order(pvalue)][1:10]

for (snp in unique(ordered_top10_snps$snp)){
    plot_PWM_heatmap_by_genotype(snpID = snp, weighting = NULL, subtitle = generate_subtitle_from_GWAS_pvalue_slope(snpID = snp, GWAS_results = ordered_top10_snps))
}
