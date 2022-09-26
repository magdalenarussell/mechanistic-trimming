OUTPUT_PATH = '/fh/fast/matsen_e/shared/tcr-gwas/mechanistic-trimming/motif_data'
PROJECT_PATH = '/home/mrussel2/mechanistic-trimming'

TCR_REPERTOIRE_DATA_parsimony = '_ignore/emerson_stats/'
TCR_REPERTOIRE_DATA_igor = '_ignore/emerson_igor_stats/'

WHOLE_NUCSEQS_parsimony = '_ignore/tcrb_processed_geneseq.tsv' 
WHOLE_NUCSEQS_igor = '_ignore/igor_imgt_genes.tsv' 
WHOLE_NUCSEQS_alpha = '_ignore/imgt_tcra_sequences.tsv' 
WHOLE_NUCSEQS_gamma = '_ignore/imgt_tcrg_sequences.tsv' 
WHOLE_NUCSEQS_delta = '_ignore/imgt_tcrd_sequences.tsv' 
WHOLE_NUCSEQS_igh = '_ignore/imgt_igh_sequences.tsv'
ALIGNED_NUCSEQS = '_ignore/combo_xcr.tsv'
CONFUSIBILITY_DATA = '_ignore/emerson_stats_ties/'

# SNP genotype file (download https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001918.v1.p1)
SNP_GDS_FILE <<- paste0(PROJECT_PATH, "/_ignore/HSCT_comb_geno_combined_v03_tcr.gds")

# File to map scanIDs to localID (data located within the `gwas_id_mapping.tsv` file available at https://doi.org/10.5281/zenodo.5719520)
#TODO change back to file without HIP ids...
# ID_MAPPING_FILE <<- paste0(PROJECT_PATH, '/_ignore/gwas_id_mapping.tsv')
ID_MAPPING_FILE <<- paste0(PROJECT_PATH, '/_ignore/temp_gwas_id_mapping.tsv')
