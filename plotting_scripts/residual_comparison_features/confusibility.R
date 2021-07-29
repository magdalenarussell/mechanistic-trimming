get_feature <- function(predicted_trims){
    files = fs::dir_ls(path = CONFUSIBILITY_DATA)
    registerDoParallel(cores=NCPU)
    together = foreach(file = files, .combine=rbind) %dopar% {
        file_data = fread(file)
        print(paste(file))
        get_confusibility_by_subject(file_data)
    }
    stopImplicitCluster()

    confus_props = together[, sum(confus_prop)/.N, by = gene]
    setnames(confus_props, 'V1', 'feature') 
    return(confus_props)
}

get_confusibility_by_subject <- function(file_data){
    repsize = nrow(file_data)
    gene_type = substring(GENE_NAME, 1, 1)
    tie_type = paste0(gene_type, '_ties')

    confus = file_data[get(tie_type) != get(GENE_NAME), .N, by = GENE_NAME]
    colnames(confus) = c('gene', 'confus')

    gene_usage = file_data[, .N, by = GENE_NAME]
    colnames(gene_usage) = c('gene', 'gene_count')

    together = merge(confus, gene_usage, all.y = TRUE)
    together[is.na(confus), confus := 0]
    together[, confus_prop := confus/gene_count]
    together$repsize = repsize
    return(together)
}

