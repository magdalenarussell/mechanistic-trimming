#I found a typo in the mclogit predict function...here is a temporary function
#to replace until it is fixed in the package
temp_predict <- function(model, newdata, se.fit = FALSE){
    formula = model$formula
    lhs = formula[[2]]
    rhs = formula[-2]
    if(length(lhs)==3)
        sets = lhs[[3]]
    else stop("no way to determine choice set ids")
    if(missing(newdata)){
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
    if(se.fit){
        variance = vcov(model)
        stopifnot(ncol(design_matrix)==ncol(variance))
    }
    set = match(set,unique(set))
    exp.eta = exp(eta)
    sum.exp.eta = rowsum(exp.eta,set)
    probabilities = exp.eta/sum.exp.eta[set]
    if(se.fit){
        wX <- probabilities*(design_matrix - rowsum(probabilities*design_matrix,set)[set,,drop=FALSE])
        se.p <- sqrt(rowSums(wX * (wX %*% variance)))
        if(is.null(na.act))
            list(fit=probabilities,se.fit=se.p)
        else
            list(fit=napredict(na.act,probabilities),
                 se.fit=napredict(na.act,se.p))
    } else {
        if(is.null(na.act)) probabilities 
        else napredict(na.act,probabilities)
    }
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
    path = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, 'temp_evaluation', type) 
    if (!dir.exists(path)){
        dir.create(path, recursive = TRUE)
    }
    return(path)
}

get_per_run_model_evaluation_file_name <- function(type){
    path = get_per_run_model_evaluation_path(type)
    name = file.path(path, paste0('model_evaluation_', MODEL_TYPE, '_', MOTIF_TYPE, '_motif_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '_bounded_', LOWER_TRIM_BOUND, '_', UPPER_TRIM_BOUND, '_', GENE_WEIGHT_TYPE, '.tsv'))
    return(name)
}

compile_result <- function(loss_list, type, parameter_count){
    loss_length = length(loss_list)
    result = data.table(motif_length_5_end = rep(LEFT_NUC_MOTIF_COUNT, loss_length), motif_length_3_end = rep(RIGHT_NUC_MOTIF_COUNT, loss_length), motif_type = rep(MOTIF_TYPE, loss_length), gene_weight_type = rep(GENE_WEIGHT_TYPE, loss_length), upper_bound = rep(UPPER_TRIM_BOUND, loss_length), lower_bound = rep(LOWER_TRIM_BOUND, loss_length), model_type = rep(MODEL_TYPE, loss_length), terminal_melting_5_end_length = rep(LEFT_SIDE_TERMINAL_MELT_LENGTH, loss_length), held_out_gene_fraction = rep(HELD_OUT_FRACTION, loss_length), sample_repetitions = rep(REPETITIONS, loss_length), model_parameter_count = parameter_count) 
    result[[type]] = loss_list
    return(result)
}
 
write_result_dt <- function(log_loss, type, parameter_count){
    file_name = get_per_run_model_evaluation_file_name(type)
    result = compile_result(log_loss, type, parameter_count)
    fwrite(result, file_name, sep = '\t')
    return(result)
}

compile_evaluation_results <- function(type){
    files = fs::dir_ls(path = get_per_run_model_evaluation_path(type))
    require(parallel)
    cluster = makeCluster(NCPU)
    files_dt = parLapply(cluster, files, function(x){
                        data.table::fread(x)})
    stopCluster(cluster)
    rbound = rbindlist(files_dt)
    return(rbound)
}

source(paste0(PROJECT_PATH, '/scripts/model_evaluation_type_functions/', TYPE, '.R'))
