
get_random_genes <- function(gene_count){
    set.seed(66)
    gene_dt = data.frame()
    for (gene in paste0('gene_', seq(1, gene_count))){
        sequence = randDNA(50) 
        gene_dt = rbind(gene_dt, c(gene, sequence))
    }
    colnames(gene_dt) = c('gene', 'whole_nucseq')
    return(gene_dt)
} 

get_gene_prob <- function(genes, gene_probs_type){
    if (gene_probs_type == 'uniform'){
        probs = rep(1, length(unique(genes$gene)))
    } else {
        probs = sample(seq(0, 1, by = 0.1), length(unique(genes$gene)), replace = TRUE)
    } 
    return(probs)
}

generate_random_sim_data <- function(subject_names, genes, gene_probs = 'uniform'){
    sim_data = data.frame()
    for (subject in subject_names){
        repsize = sample(100:1000, 1, replace = TRUE)
        probs = get_gene_prob(genes, gene_probs)
        for (tcr in seq(repsize)){
            gene_sample = sample(1:nrow(genes), 1, replace = TRUE, prob = probs)
            gene = genes[gene_sample,]
            tcr_obs = c(subject = subject, gene)
            sim_data = rbind(sim_data, tcr_obs)
        }
    }
    return(as.data.table(sim_data))
}

create_position_base_variables <- function(dataframe){
    positions = get_positions()
    for (position in positions){
        for (base in c('A', 'T', 'C', 'G')){
            dataframe[get(position) == base, paste0(position, base) := 1]
            dataframe[get(position) != base, paste0(position, base) := 0]
        }
    }
    return(dataframe)
}

get_position_count_matrix <- function(data){
    total_positions = LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT 
    positions = get_positions()
    pcm = matrix(0, ncol = total_positions, nrow = 4)
    rownames(pcm) = c('A', 'C', 'T', 'G')
    colnames(pcm) = positions 

    for (position in positions){
        for (base in c('A', 'T', 'C', 'G')){
            counts = sum(data[[paste0(position, base)]])
            pcm[base, position] = counts
        }
    }
    return(pcm)
}

get_position_probability_matrix <- function(pcm){
    ppm = pcm/colSums(pcm)
    return(ppm)
}

get_background_frequencies_by_gene <- function(data){
    sequences = unique(data$whole_nucseq)
    counts = matrix(0, ncol = 1, nrow = 4)
    rownames(counts) = c('A', 'C', 'T', 'G')
    for (base in c('A', 'T', 'C', 'G')){
        counts[base, 1] = 1/(sum(str_count(sequences, base))/sum(nchar(sequences)))
    }
    return(counts)
}

get_pwm <- function(data){
    require(matrixcalc)
    positions = get_positions()
    total_positions = LEFT_NUC_MOTIF_COUNT + RIGHT_NUC_MOTIF_COUNT 
    pcm = get_position_count_matrix(data)
    ppm = get_position_probability_matrix(pcm)
    background = get_background_frequencies_by_gene(data)
    background_matrix = matrix(rep(background, total_positions), 4)
    rownames(background_matrix) = c('A', 'C', 'T', 'G')
    colnames(background_matrix) = positions 
    pwm = log(hadamard.prod(ppm, background_matrix))
    return(pwm)
}

get_gene_trim_dataframe <- function(gene_count, trim_lengths){
    genes_trim = data.table()
    genes = get_random_genes(gene_count)

    for (trim in trim_lengths){
        genes$trim_length = trim
        genes$count = 100
        genes_trim = rbind(genes_trim, genes)
    }
    genes_trim$motif = get_motif_context_unobserved(genes_trim$whole_nucseq, genes_trim$trim_length)
    final = split_motif_column_by_motif_position(genes_trim)

    return(final)
}

generate_random_PWM <- function(gene_count){
    trims = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)

    genes_trim = get_gene_trim_dataframe(gene_count, trims)
    
    final = create_position_base_variables(genes_trim)

    PWM = get_pwm(final)

    return(PWM)
}

get_pwm_weight_single_trim <- function(PWM, gene_trim_df_row){
    sum = 0
    positions = get_positions()
    for (position in positions){
        sum = sum + PWM[gene_trim_df_row[[position]], position]
    }
    return(exp(sum))
}

get_probability_from_pwm <- function(PWM, gene_trim_df){
    single_pwm_weights = c()
    for (row in seq(nrow(gene_trim_df))){
        weight = get_pwm_weight_single_trim(PWM, gene_trim_df[row,])
        single_pwm_weights = c(single_pwm_weights, weight)
    }

    gene_trim_df$single_position_weight = single_pwm_weights
    gene_trim_df[, all_position_weight_sum := sum(single_position_weight), by = gene]
    gene_trim_df$p_trim_given_gene_pwm = gene_trim_df$single_position_weight/gene_trim_df$all_position_weight_sum 

    return(gene_trim_df[, c('gene', 'whole_nucseq', 'trim_length', 'motif', 'p_trim_given_gene_pwm')])
}


get_trimming_amounts_using_PWM <- function(genes, sim_data, PWM){
    trims = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND)

    gene_trim_df = get_gene_trim_dataframe(length(unique(genes$gene)), trims)    
    probabilities = get_probability_from_pwm(PWM, gene_trim_df)

    complete_data = data.table()
    for (row in seq(nrow(sim_data))){
        temp_row = sim_data[row,]
        temp_gene = temp_row$gene
        probs = probabilities[gene == temp_gene]
        trim = sample(probs$trim_length, 1, prob = probs$p_trim_given_gene_pwm)
        temp_row$trim_length = trim 
        complete_data = rbind(complete_data, temp_row)
    }

    return(complete_data)
}

get_file_path <- function(gene_count, gene_probs_type){
    dir = paste0(gene_probs_type, '_', gene_count, '_genes')
    path = file.path(SIM_OUTPUT_LOC, dir)
    dir.create(path, recursive = TRUE)
    return(path)
}

get_file_name <- function(subject_names, gene_count, gene_probs_type){
    subj_concat = paste(subject_names, collapse = '')
    subject_count = length(subject_names)
    name = paste0(subj_concat, '_simulation_data_', subject_count, '_subject_', gene_count, '_genes_with_', gene_probs_type, '_prob.tsv')
    path = get_file_path(gene_count, gene_probs_type)
    final_loc = file.path(path, name)
    return(final_loc)
}
 

get_complete_sim_data <- function(subject_names, gene_count, PWM, gene_probs = 'uniform', gene_probs_type = 'uniform'){
    genes = get_random_genes(gene_count)
    sim_data = generate_random_sim_data(subject_names, genes, gene_probs)
    sim_data = get_trimming_amounts_using_PWM(genes, sim_data, PWM)
    sim_data$motif = get_motif_context_unobserved(sim_data$whole_nucseq, sim_data$trim)
    sim_data$gene_type = 'simulation'
    sim_data$observed = TRUE
    file_name = get_file_name(subject_names, gene_count, gene_probs_type)
    fwrite(sim_data, file = file_name, sep = '\t')

    return(sim_data)
}
