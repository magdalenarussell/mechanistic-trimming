library(foreach)
library(doParallel)
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

# compile PFM for whole population
compile_PFMs(group_localID_list = get_subject_partition_by_SNP(), group_name = get_group_name())

# compile PFM for each significant artemis SNP by genotype
significant_artemis_v_trim = compile_snps_from_GWAS(phenotype_list = c('v_trim'), gene = 'artemis') 
artemis_snp_genotypes = get_genotypes_from_GWAS_results(significant_artemis_v_trim)

for (snpID in paste(unique(significant_artemis_v_trim$snp))){
    for (genotype in c(0, 1, 2)){  

        compile_PFMs(group_localID_list = get_subject_partition_by_SNP(snpID, genotypes_df = artemis_snp_genotypes, genotype), group_name = get_group_name(snpID, genotype))
    }
}
