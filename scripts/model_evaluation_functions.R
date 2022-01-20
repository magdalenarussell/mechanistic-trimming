generate_hold_out_sample <- function(motif_data, sample_size){
    stopifnot(sample_size < length(unique(motif_data$gene)))
    # sample size is the number of gene, subject combos
    sample_genes = sample(unique(motif_data$gene), sample_size, replace = FALSE)
    sample_data = data.table()
    motif_data_subset = motif_data
    for (sample in sample_genes){
        temp = motif_data[gene == sample]
        sample_data = rbind(sample_data, temp)
        cols = colnames(sample_data)
        motif_data_subset = motif_data_subset[!sample_data, on = cols]
    }
    # updating the total_tcr, p_gene, gene_weight_type, and weighted_observation variables for the newly sampled datasets
    source(paste0('scripts/sampling_procedure_functions/', GENE_WEIGHT_TYPE, '.R'), local = TRUE)
    motif_data_subset = calculate_subject_gene_weight(motif_data_subset)
    source(paste0('scripts/sampling_procedure_functions/p_gene_given_subject.R'), local = TRUE)
    sample_data = calculate_subject_gene_weight(sample_data)
    return(list(sample = sample_data, motif_data_subset = motif_data_subset))
}

get_hold_out_sample_probability <- function(sample_size, motif_data){
    total_genes = length(unique(motif_data$gene))
    prob = (1/total_genes)^sample_size
    return(prob)
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

calculate_cond_expected_log_loss <- function(model, sample_data){
    #TODO switch this function to mclogit.predict (if typo is fixed..."contasts")
    sample_data = process_data_for_model_fit(sample_data)
    sample_data$prediction = temp_predict(model, newdata = sample_data)
    sample_data[, log_prediction := log(prediction)]
    sample_data[, weighted_log_prediction := log_prediction * weighted_observation] 
    log_loss = -sum(sample_data$weighted_log_prediction)
    return(log_loss)
}

evaluate_cond_log_loss <- function(motif_data, held_out_fraction, repetitions, write_intermediate_loss = FALSE) {
    set.seed(66)
    gene_count = length(unique(motif_data$gene))
    held_out_gene_count = round(held_out_fraction*gene_count)
    log_loss_vector = c()
    sample_prob_vector = c()
    for (rep in 1:repetitions){
        # Generate a held out sample and motif data subset
        sample_data = generate_hold_out_sample(motif_data, sample_size = held_out_gene_count) 
        motif_data_subset = sample_data$motif_data_subset
        sample = sample_data$sample

        # Fit model to the motif_data_subset
        model = fit_model(motif_data_subset)

        # Compute conditional logistic loss value for held out sample using model
        log_loss = calculate_cond_expected_log_loss(model, sample)
        log_loss_vector = c(log_loss_vector, log_loss)

        # Compute probability of held out sample
        prob = get_hold_out_sample_probability(held_out_gene_count, motif_data)
        sample_prob_vector = c(sample_prob_vector, prob)
        if (isTRUE(write_intermediate_loss)){
            write_result_dt(log_loss, type = 'log_loss', held_out_gene_fraction = held_out_fraction, repetitions = paste0('rep_', rep), intermediate = TRUE) 
        }
    }

    #TODO: should this be just multiplied by the sample_prob_vector or be the mean? 
    # expected_log_loss = sum(sample_prob_vector * log_loss_vector)
    expected_log_loss = sum((1/repetitions) * log_loss_vector)

    return(expected_log_loss)
}

get_model_evaluation_file_name <- function(type, intermediate){
    stopifnot(type %in% c('log_loss'))
    path = file.path(OUTPUT_PATH, ANNOTATION_TYPE, TRIM_TYPE, PRODUCTIVITY) 
    if (isTRUE(intermediate)){
        name = file.path(path, paste0('intermediate_model_evaluation_expected_', type, '.tsv'))
    } else {
        name = file.path(path, paste0('model_evaluation_expected_', type, '.tsv'))
    }
    return(name)
}

write_result_dt <- function(log_loss, type, held_out_gene_fraction, repetitions, intermediate = FALSE){
    file_name = get_model_evaluation_file_name(type, intermediate)
    result = data.table(motif_length_5_end = LEFT_NUC_MOTIF_COUNT, motif_length_3_end = RIGHT_NUC_MOTIF_COUNT, motif_type = MOTIF_TYPE, gene_weight_type = GENE_WEIGHT_TYPE, upper_bound = UPPER_TRIM_BOUND, lower_bound = LOWER_TRIM_BOUND, model_type = MODEL_TYPE, terminal_melting_5_end_length = LEFT_SIDE_TERMINAL_MELT_LENGTH, held_out_gene_fraction = held_out_gene_fraction, sample_repetitions = repetitions) 
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
