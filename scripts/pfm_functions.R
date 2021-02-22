source(paste0('scripts/motif_class_functions/', MOTIF_TYPE, '.R'))
source(paste0('scripts/pfm_class_functions/', PFM_TYPE, '.R'))
source('scripts/GWAS_scripts/regression_parameters.R')
source('scripts/GWAS_scripts/file_paths.R')

extract_subject_ID <- function(tcr_repertoire_file_path){
    file_name = str_split(tcr_repertoire_file_path, "/")[[1]][3]
    file_root_name = str_split(file_name, ".tsv")[[1]][1]
    localID = str_split(file_root_name, "_")[[1]][3]
    return(localID)
}

# The following two functions were designed to calculate PFMs by trim_length
# compute_PFMs_for_all_trim_lengths <- function(data, file_path){
    # for (trim_length in unique(data[[TRIM_TYPE]])){
    #     if (trim_length < 1){
    #         next
    #     }
    #     trim_length_PFM = as.data.frame(create_PFM_by_trim_length(trim_length, data))
    #     colnames(trim_length_PFM) = c(paste0(seq(-LEFT_NUC_MOTIF_COUNT,-1)), paste0(seq(1, RIGHT_NUC_MOTIF_COUNT)))
    #     write.table(trim_length_PFM, file = paste0(file.path(file_path, paste0('PFM_trim_', as.character(trim_length), '.tsv'))), sep = '\t') 
    # }
# }

# create_PFMs_for_subject <- function(file_path){
    # temp_data = fread(file_path)
    # # TODO, for now, don't remove missing d gene instances
    # # temp_data = temp_data[d_gene != '-']
    # subject_id = extract_subject_ID(file_path)
    
    # output_path = file.path('_ignore', 'pfms', paste0(PFM_TYPE, '_', MOTIF_TYPE), paste0(subject_id))
    # dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
    # whole_nucseq = fread('_ignore/tcrb_processed_geneseq.tsv')
    # temp_data = merge(temp_data, whole_nucseq, by.x = GENE_NAME, by.y = 'gene_names')

    # compute_PFMs_for_all_trim_lengths(temp_data, output_path)
# }

compute_PFMs_for_all_genes <- function(data, file_path){
    for (gene in unique(data[[GENE_NAME]])){
        data_subset = data[get(GENE_NAME)==gene]
        gene_PFM = as.data.frame(matrix(0, ncol = LEFT_NUC_MOTIF_COUNT+RIGHT_NUC_MOTIF_COUNT, nrow = 5))
        for (trim_length in unique(data_subset[[TRIM_TYPE]])){
            if (trim_length < 1){
                trim_length_PFM = as.data.frame(matrix(0, ncol = LEFT_NUC_MOTIF_COUNT+RIGHT_NUC_MOTIF_COUNT, nrow = 5))
                rownames(trim_length_PFM) = c('A', 'C', 'G', 'T', 'other')
            } else {
                trim_length_PFM = as.data.frame(create_PFM_by_trim_length(trim_length, data_subset))
            }
            colnames(trim_length_PFM) = c(paste0(seq(LEFT_NUC_MOTIF_COUNT,1)), paste0(seq(-1, -RIGHT_NUC_MOTIF_COUNT)))
            gene_PFM = trim_length_PFM + gene_PFM
        }
        if (length(str_split(gene, '/')[[1]]) > 1){
            gene = str_replace(gene, '/', '_')
        }
        write.table(gene_PFM, file = paste0(file.path(file_path, paste0('PFM_gene_', gene, '.tsv'))), sep = '\t') 
    }
}

create_PFMs_for_subject <- function(file_path){
    temp_data = fread(file_path)
    # TODO, for now, don't remove missing d gene instances
    # temp_data = temp_data[d_gene != '-']
    subject_id = extract_subject_ID(file_path)
    
    output_path = file.path('_ignore', 'pfms', paste0(PFM_TYPE, '_', MOTIF_TYPE), paste0(subject_id))
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
    whole_nucseq = fread('_ignore/tcrb_processed_geneseq.tsv')
    temp_data = merge(temp_data, whole_nucseq, by.x = GENE_NAME, by.y = 'gene_names')

    compute_PFMs_for_all_genes(temp_data, output_path)
}



create_all_PFMs <- function(directory){
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    foreach(file = files) %dopar% {
        create_PFMs_for_subject(file)
        print(paste0(file))
    }
    stopImplicitCluster()
}

# This function compiled PFMs calculated by trim_length
# compile_PFMs <- function(group_localID_list, group_name, trim_lengths_list = 'ALL'){
#     compiled_PFM = matrix(0, ncol = (RIGHT_NUC_MOTIF_COUNT+LEFT_NUC_MOTIF_COUNT), nrow = 5)
#     trim_length_file_names = paste0('PFM_trim_', trim_lengths_list, '.tsv')
#     for (subject in group_localID_list){
#         path = file.path('_ignore', 'pfms', paste0(PFM_TYPE, '_', MOTIF_TYPE), paste0(subject))
#         files = fs::dir_ls(path = path)
#         if (trim_lengths_list != 'ALL'){
#             files = files[files %in% paste0(path, '/', trim_length_file_names)]
#         }
#         for (file in files){
#             temp_matrix = as.matrix(fread(file), rownames = 1)
#             compiled_PFM = compiled_PFM + temp_matrix
#         }
#     }

#     output_path = file.path('_ignore', 'pfms', 'complete') 
#     dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
#     write.table(as.data.frame(compiled_PFM), file = file.path(output_path, paste0(group_name, '_', PFM_TYPE, '_', MOTIF_TYPE, '_', paste(trim_lengths_list, collapse = ''), '.tsv')), sep = '\t')
# }

compile_PFMs <- function(group_localID_list, group_name, weighting = NULL){
    compiled_PFM = matrix(0, ncol = (RIGHT_NUC_MOTIF_COUNT+LEFT_NUC_MOTIF_COUNT), nrow = 5)
    for (subject in group_localID_list){
        path = file.path('_ignore', 'pfms', paste0(PFM_TYPE, '_', MOTIF_TYPE), paste0(subject))
        files = fs::dir_ls(path = path)
        for (file in files){
            if (!is.null(weighting)){
                #TODO add weighting scheme
            }
            temp_matrix = as.matrix(fread(file), rownames = 1)
            compiled_PFM = compiled_PFM + temp_matrix
        }
    }

    output_path = file.path('_ignore', 'pfms', 'complete') 
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
    write.table(as.data.frame(compiled_PFM), file = file.path(output_path, paste0(group_name, '_', PFM_TYPE, '_', MOTIF_TYPE, '_', weighting, '.tsv')), sep = '\t')
}


find_base_background_frequencies <- function(){
    #TODO not totally sure what to do here, but I think just find the
    #frequencies across ALL possible v-genes for example (but maybe I should be
    #doing this on a per-subject basis...)
    #OR compute an A/T and G/C background proportion...since pnucs deal with
    #either or...
    whole_nucseq = fread('_ignore/tcrb_processed_geneseq.tsv')
    whole_nucseq$gene_type = paste0(tolower(substr(whole_nucseq$gene_names, 4, 4)), '_gene')
    whole_nucseq_by_gene = whole_nucseq[gene_type == GENE_NAME]
    sequence_composition = letterFrequency(DNAStringSet(whole_nucseq_by_gene$sequences), letters = c('T', 'G', 'C', 'A'))
    background_frequency = colSums(sequence_composition)/sum(sequence_composition)
    return(background_frequency)
}

calculate_PPM_from_PFM <- function(PFM){
    PFM_filtered = PFM[c('T', 'G', 'C', 'A'),]
    PPM = t(t(PFM_filtered)/colSums(PFM_filtered))
    return(PPM)
}

calculate_PWM_from_PFM <- function(PFM, weighting){
    background_freqs = find_base_background_frequencies()
    background_matrix = matrix(background_freqs, length(background_freqs), (RIGHT_NUC_MOTIF_COUNT+LEFT_NUC_MOTIF_COUNT))

    PPM = calculate_PPM_from_PFM(PFM)
    rownames(background_matrix) = names(background_freqs)
    stopifnot(identical(rownames(background_matrix), rownames(PPM)))

    PWM = log(PPM/background_matrix)
    return(PWM)
}

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

