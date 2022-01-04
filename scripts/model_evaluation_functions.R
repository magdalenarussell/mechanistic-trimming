generate_hold_out_sample <- function(motif_data, sample_size){
    stopifnot(sample_size < length(unique(motif_data$subject))*length(unique(motif_data$gene)))
    # sample size is the number of gene, subject combos
    set.seed(66)
    sample_genes = sample(unique(motif_data$gene), sample_size, replace = TRUE)
    sample_data = data.table()
    motif_data_subset = motif_data
    for (sample in sample_genes){
        sample_subject = sample(unique(motif_data$subject), 1)
        temp = motif_data[subject == sample_subject & gene == sample]
        sample_data = rbind(sample_data, temp)
        cols = colnames(sample_data)
        motif_data_subset = motif_data_subset[!sample_data, on = cols]
    }
    return(list(sample = sample_data, motif_data_subset = motif_data_subset))
}

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

calculate_cond_log_loss <- function(model, sample_data){
    #TODO switch this function to mclogit.predict (if typo is fixed..."contasts")
    sample_data = process_data_for_model_fit(sample_data)
    sample_data$prediction = temp_predict(model, newdata = sample_data)
    sample_data$log_prediction = log(sample_data$prediction)
    sample_data$weighted_log_prediction = sample_data$log_prediction * sample_data$count
    log_loss = -sum(sample_data$weighted_log_prediction)
    return(log_loss)
}

get_model_evaluation_file_name <- function(type){
    stopifnot(type %in% c('log_loss', 'per_gene', 'per_gene_per_trim'))
    name = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY, paste0('model_evaluation_', type, '.tsv'))
    return(name)
}

write_result_dt <- function(log_loss, type){
    file_name = get_model_evaluation_file_name(type)
    result = data.table(motif_length_5_end = LEFT_NUC_MOTIF_COUNT, motif_length_3_end = RIGHT_NUC_MOTIF_COUNT, motif_type = MOTIF_TYPE, gene_weight_type = GENE_WEIGHT_TYPE, upper_bound = UPPER_TRIM_BOUND, lower_bound = LOWER_TRIM_BOUND, model_type = MODEL_TYPE, terminal_melting_5_end_length = LEFT_SIDE_TERMINAL_MELT_LENGTH) 
    result[[type]] = log_loss
    
    if (file.exists(file_name)){
        results = fread(file_name)
        together = rbind(results, result)
        fwrite(together, file_name, sep = '\t')
    } else {
        fwrite(result, file_name, sep = '\t')
    }
    return(result)
}

evaluate_model_per_gene <- function(data, type){
    if (type == 'per_gene'){
        rmse = calculate_rmse_by_gene(data)
    } else if (type == 'per_gene_per_trim'){
        rmse = calculate_rmse(data)
    }

    mean_abs_resid = mean(rmse$rmse)
    return(mean_abs_resid)
}
