get_model_formula <- function(){
    formula = formula(paste0('cbind(weighted_observation, interaction(gene, subject)) ~ as.factor(trim_length)'))
    return(formula)
}


