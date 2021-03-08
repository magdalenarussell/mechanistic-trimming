source(paste0('scripts/motif_class_functions/', MOTIF_TYPE, '.R'))
source(paste0('scripts/pfm_class_functions/', PFM_TYPE, '.R'))

extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][3]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_")[[1]][3]
    return(localID)
}

get_subject_motif_output_location <- function(){
    output_location = file.path(OUTPUT_PATH, paste0(MOTIF_TYPE, '_', PFM_TYPE))
    return(output_location)
}

compile_motifs_for_subject <- function(file_path){
    temp_data = fread(file_path)
    # TODO, for now, don't remove missing d gene instances
    # temp_data = temp_data[d_gene != '-']
    subject_id = extract_subject_ID(file_path)
    
    output_location = get_subject_motif_output_location() 
    dir.create(output_location, recursive = TRUE, showWarnings = FALSE)
    whole_nucseq = fread('_ignore/tcrb_processed_geneseq.tsv')
    temp_data = merge(temp_data, whole_nucseq, by.x = GENE_NAME, by.y = 'gene_names')

    motif_data = get_motifs(temp_data, subject_id)
    fwrite(motif_data, file = file.path(output_location, paste0(subject_id, '.tsv')), sep = '\t')
}


compile_all_motifs <- function(directory){
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    foreach(file = files) %dopar% {
        compile_motifs_for_subject(file)
        print(paste0(file))
    }
    stopImplicitCluster()
}

calculate_pdel_seq_and_gene <- function(motif_data){
    motif_data = as.data.table(motif_data)
    motif_data[,count_subject_gene_seq := .N, by = .(gene, gene_type, subject, trim_length)]
    motif_data[,count_subject := .N, by = .(subject)]
    motif_data[observed == TRUE,pdel_seq_and_gene := count_subject_gene_seq/count_subject]
    motif_data[observed == FALSE, pdel_seq_and_gene := 0]
    return(motif_data)
}

calculate_pdel_seq_and_gene2 <- function(motif_data){
    motif_data = as.data.table(motif_data)
    motif_data = motif_data[observed == TRUE]
    motif_data[,count_subject_gene_trim_length := .N, by = .(gene, gene_type, subject, trim_length)]
    motif_data[,count_subject := .N, by = .(subject)]
    motif_data[, count_subject_gene := .N, by = .(gene, subject)]
    motif_data[, p_trim_and_gene := count_subject_gene_trim_length/count_subject]
    motif_data[, p_gene := count_subject_gene/count_subject]
    motif_data[, p_trim_given_gene := p_trim_and_gene/p_gene]
    return(unique(motif_data))
}


split_motif_column_by_motif_position <- function(motif_dataframe){
    motif_dataframe = motif_dataframe %>% separate(motif, c(paste0('pos5_', c(seq(LEFT_NUC_MOTIF_COUNT, 1))), paste0('pos3_', c(seq(1, RIGHT_NUC_MOTIF_COUNT)))), sep = seq(1, LEFT_NUC_MOTIF_COUNT+RIGHT_NUC_MOTIF_COUNT-1))
    return(motif_dataframe)
}

prepare_data_for_regression <- function(motif_data){
    motif_data = calculate_pdel_seq_and_gene2(motif_data)
    motif_data = split_motif_column_by_motif_position(motif_data)
    motif_data[motif_data == '-'] <- NA
    motif_data = unique(motif_data)
    return(motif_data)
}

concatenate_motifs_all_subjects <- function(directory = get_subject_motif_output_location()){
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    together = foreach(file = files, .combine=rbind) %dopar% {
        file_data = fread(file)
        print(paste(file))
        prepare_data_for_regression(file_data)
    }
    stopImplicitCluster()
    return(together)
}
    
source(paste0('scripts/regression_class_functions/', REGRESSION_TYPE, '_model.R'))

compile_snps_from_GWAS <- function(phenotype_list, gene, pvalue_cutoff = 0.05/35481497){
    compiled_data = data.table()
    for (phenotype in phenotype_list){
        parameters = set_regression_parameters(phenotype)
        file_name = get_compiled_file_name_with_arguments(parameters[['phenotype']], parameters[['condensing_variable']], parameters[['infer_missing_d_gene']], parameters[['pca_count']])
        compiled_data = rbind(compiled_data, fread(file_name))
    }
    
    genes = c('artemis', 'mhc', 'dntt', 'rag', 'tcrb', 'tcra')
    chr = c(10, 6, 10, 11, 7, 14)
    pos1 = c(14939358, 25912984, 98064085, 36510709, 141998851, 22090057)
    pos2 = c(14996431, 33290793, 98098321, 36593156, 142510972, 23021075)

    gene_annotations = data.table(genes = genes, chr = chr, pos1 = pos1, pos2 = pos2)
    gene_coords = gene_annotations[genes == gene]
    gene_subset_compiled_data = compiled_data[hg19_pos < (gene_coords$pos2 + 200000) & hg19_pos > (gene_coords$pos1 - 200000) & chr == gene_coords$chr & pvalue < pvalue_cutoff]
    return(gene_subset_compiled_data)
}

map_scanID_to_localID <- function(scanIDs_to_convert){
    ID_map_file = fread(ID_MAPPING_FILE)
    converted_IDs = plyr::mapvalues(scanIDs_to_convert,
                                    ID_map_file$scanID,
                                    ID_map_file$localID)
    return(converted_IDs)
}

get_genotypes_by_snpID <- function(snpID_list){
    snp_gds_file = snpgdsOpen(SNP_GDS_FILE)
    genotypes = snpgdsGetGeno(snp_gds_file, snp.id = snpID_list, with.id = TRUE)
    genotypes_df = as.data.frame(genotypes$genotype)
    genotypes_df$scanID = genotypes$sample.id
    colnames(genotypes_df) = c(paste0(genotypes$snp.id), 'scanID')
    snpgdsClose(snp_gds_file)    
    
    genotypes_df$localID = map_scanID_to_localID(genotypes_df$scanID)
    return(genotypes_df)
}

get_genotypes_from_GWAS_results <- function(GWAS_results){
    snps = unique(GWAS_results$snp)
    genotypes_df = get_genotypes_by_snpID(snps)
    return(genotypes_df)
}

get_subject_partition_by_SNP <- function(snpID = NA, genotypes_df = NA, genotype = NA){
    if (is.na(genotype)){
        subjects = list.dirs(path = paste0("_ignore/pfms/", PFM_TYPE, '_', MOTIF_TYPE), full.names= FALSE, recursive = FALSE)
    } else {
        stopifnot(genotype %in% c(0,1,2))
        stopifnot(!is.na(snpID))
        genotypes_dt = as.data.table(genotypes_df)
        subjects = genotypes_dt[get(snpID) == genotype]$localID
    }
    return(subjects)
}

get_group_name <- function(snpID = NA, genotype = NA){
    if (is.na(snpID) & is.na(genotype)){
        group_name = 'ALL'
    } else {
        stopifnot(genotype %in% c(0,1,2))
        stopifnot(!is.na(snpID))
        group_name = paste0(snpID, '_genotype_', genotype)
    }
    return(group_name)
}

