#I found a typo in the mclogit predict function...here is a temporary function
#to replace until it is fixed in the package
temp_predict <- function(model, newdata, se.fit = FALSE){
    formula = model$formula
    lhs = formula[[2]]
    rhs = formula[-2]
    if (length(lhs)==3) {
        sets = lhs[[3]]
    } else {
        stop("no way to determine choice set ids")
    }
    if (missing(newdata)){
        model_data =  model.frame(formula,data=model$data)
        set = model_data[[1]][,2]
        na.act = model$na.action
    } else{
        lhs = lhs[[3]]
        formula[[2]] = lhs
        model_data = model.frame(formula,data=newdata)
        set = model_data[[1]]
        na.act = attr(model_data,"na.action")
    }
    design_matrix = model.matrix(rhs,model_data,
                      contrasts.arg=model$contrasts,
                      xlev=model$xlevels
                      )
    coefs = coef(model)
    design_matrix = design_matrix[,names(coefs), drop=FALSE]
        
    eta = c(design_matrix %*% coefs)
    if (se.fit){
        variance = vcov(model)
        stopifnot(ncol(design_matrix)==ncol(variance))
    }
    set = match(set,unique(set))
    exp.eta = exp(eta)
    sum.exp.eta = rowsum(exp.eta,set)
    probabilities = exp.eta/sum.exp.eta[set]
    if (se.fit){
        wX <- probabilities*(design_matrix - rowsum(probabilities*design_matrix,set)[set,,drop=FALSE])
        se.p <- sqrt(rowSums(wX * (wX %*% variance)))
        if (is.null(na.act))
            list(fit=probabilities,se.fit=se.p)
        else
            list(fit=napredict(na.act,probabilities),
                 se.fit=napredict(na.act,se.p))
    } else {
        if (is.null(na.act)) probabilities 
        else napredict(na.act,probabilities)
    }
}

aggregate_validation_data <- function(directory, trim_type = TRIM_TYPE, gene_type = GENE_NAME){
    files = fs::dir_ls(path = directory)
    registerDoParallel(cores=NCPU)
    together = foreach(file = files, .combine=rbind) %dopar% {
        print(paste(file))
        temp = compile_data_for_subject(file, write = FALSE, gene_type = gene_type, trim_type = trim_type)
        if (nrow(temp) > 0) {
            temp
        }
    }
    
    trims = get_trim_order(trim_type)
    genes = get_gene_order(gene_type)

    if (MODEL_TYPE %like% 'dna_shape') {
        together = convert_data_to_motifs(together, left_window_size = LEFT_NUC_MOTIF_COUNT + 2, right_window_size = RIGHT_NUC_MOTIF_COUNT + 2)
        processed_motif_data = process_data_for_model_fit(together, gene_type = gene_type, trim_type = trim_type)
        motif_data = convert_data_to_motifs(processed_motif_data)
    } else {
        together = convert_data_to_motifs(together)
        motif_data = process_data_for_model_fit(together, gene_type = gene_type, trim_type = trim_type)
    }

    cols = colnames(motif_data)[!(colnames(motif_data) %like% 'left_nucs')]
    cols = cols[!(cols %like% 'right_nucs')]

    #remove genes with ambiguous bases
    motvar = paste0(trims, '_motif')
    genevar = paste0(genes, '_group')
    amb = c()
    for (i in seq(length(genevar))){
        temp = unique(motif_data[get(motvar[i]) %like% 'N'][[paste0(genevar[i])]])
        amb = c(amb, temp)
        motif_data = motif_data[!(get(paste0(genevar[i])) %in% amb)]
    }

    motif_data = motif_data[, ..cols]
    together_pos = split_motif_column_by_motif_position(motif_data, trim_type) 
    weighted_together = calculate_subject_gene_weight(together_pos, gene_type = gene_type, trim_type = trim_type)
    stopImplicitCluster()
    return(weighted_together)
}

calculate_cond_expected_log_loss <- function(model, sample_data){
    #TODO switch this function to mclogit.predict (if typo is fixed..."contasts")
    if (MODEL_TYPE != 'null') {
        sample_data$prediction = temp_predict(model, newdata = sample_data)
    } else {
        sample_data$prediction = exp(1)/sum(rep(exp(1), UPPER_TRIM_BOUND -LOWER_TRIM_BOUND + 1))
    }
    sample_data[, log_prediction := log(prediction)]
    sample_data[, weighted_log_prediction := log_prediction * weighted_observation] 
    log_loss = -sum(sample_data$weighted_log_prediction)
    return(log_loss)
}

get_per_run_model_evaluation_path <- function(type){
    path = file.path(MOD_OUTPUT_PATH, ANNOTATION_TYPE, PARAM_GROUP, 'temp_evaluation', LOSS_GENE_WEIGHT, type) 
    if (!dir.exists(path)){
        dir.create(path, recursive = TRUE)
    }
    if (type == 'expected_log_loss'){
        path2 = paste0(path, '/raw')
        if (!dir.exists(path2)){
            dir.create(path2, recursive = TRUE)
        }
    }
    return(path)
}

get_per_run_model_evaluation_file_name <- function(type){
    path = get_per_run_model_evaluation_path(type)
    name = file.path(path, paste0('model_evaluation_', MODEL_TYPE, '_', MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND, '_', GENE_WEIGHT_TYPE, '.tsv'))
    return(name)
}

compile_result <- function(loss_list, type, parameter_count, held_out_genes = NA, held_out_clusters = NA, validation_gene_weighting = NA){
    loss_length = length(loss_list)
    result = data.table(motif_length_5_end = rep(LEFT_NUC_MOTIF_COUNT, loss_length), motif_length_3_end = rep(RIGHT_NUC_MOTIF_COUNT, loss_length), motif_type = rep(MOTIF_TYPE, loss_length), gene_weight_type = rep(GENE_WEIGHT_TYPE, loss_length), upper_bound = rep(UPPER_TRIM_BOUND, loss_length), lower_bound = rep(LOWER_TRIM_BOUND, loss_length), model_type = rep(MODEL_TYPE, loss_length), terminal_melting_5_end_length = rep(LEFT_SIDE_TERMINAL_MELT_LENGTH, loss_length), held_out_gene_fraction = rep(HELD_OUT_FRACTION, loss_length), sample_repetitions = rep(REPETITIONS, loss_length), model_parameter_count = parameter_count, held_out_genes = held_out_genes, held_out_clusters = held_out_clusters, loss_gene_weighting = validation_gene_weighting) 
    result[[unique(type)]] = loss_list
    return(result)
}
 
write_result_dt <- function(log_loss, type, parameter_count, held_out_clusters, held_out_genes, validation_gene_weighting = NA){
    file_name = get_per_run_model_evaluation_file_name(unique(type))
    result = compile_result(log_loss, type, parameter_count, held_out_genes, held_out_clusters, validation_gene_weighting)
    fwrite(result, file_name, sep = '\t')
    return(result)
}

compile_evaluation_results <- function(type){
    files = fs::dir_ls(path = get_per_run_model_evaluation_path(type), glob = '*.tsv')
    require(parallel)
    cluster = makeCluster(NCPU)
    files_dt = parLapply(cluster, files, function(x){
                        data.table::fread(x)})
    stopCluster(cluster)
    rbound = rbindlist(files_dt, fill = TRUE)
    rbound$loss_type = type
    return(rbound)
}

source(paste0(MOD_PROJECT_PATH, '/scripts/model_evaluation_type_functions/', TYPE, '.R'))
