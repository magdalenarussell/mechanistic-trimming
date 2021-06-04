library(RBioinf)

get_random_gene_sequence <- function(motif_length, min_trim_length, max_trim_length){
    trim_count = max_trim_length - min_trim_length
    seq_length = trim_count + motif_length
    seq = randDNA(seq_length)
    return(seq)
}

get_motif_from_gene_sequence_trim <- function(gene_sequence, motif_length, trim_length){
    length = nchar(gene_sequence)
    motif = substring(gene_sequence, (length - motif_length + (motif_length/2 - trim_length) + 1), length + (motif_length/2 - trim_length))
    return(motif)
}

separate_motif_column_by_position <- function(dataframe){
    split = dataframe %>% 
        separate(motif, paste0('pos', seq(1, 8)), sep = seq(1, 8-1))
    return(split)
}

create_position_base_variables <- function(dataframe){
    for (position in paste0('pos', seq(1, 8))){
        for (base in c('A', 'T', 'C', 'G')){
            dataframe[get(position) == base, paste0(position, base) := 1]
            dataframe[get(position) != base, paste0(position, base) := 0]
        }
    }
    return(dataframe)
}

get_probability_from_pwm <- function(pwm, dataframe){
    log_prob = c()
    for (row in seq(nrow(dataframe))){
        sum = 0
        for (position in colnames(pwm)){
            sum = sum + pwm[dataframe[row][[position]], position]
        }
        log_prob = c(log_prob, sum)
    }
    dataframe$pwm_log_prob = log_prob
    dataframe$pwm_prob = exp(dataframe$pwm_log_prob)
    return(dataframe)
}

get_position_count_matrix <- function(data){
    pcm = matrix(0, ncol = 8, nrow = 4)
    rownames(pcm) = c('A', 'C', 'T', 'G')
    colnames(pcm) = paste0('pos', seq(1, 8))

    for (position in paste0('pos', seq(1, 8))){
        for (base in c('A', 'T', 'C', 'G')){
            counts = data$count * data[[paste0(position, base)]]
            total = sum(counts)
            pcm[base, position] = total
        }
    }
    return(pcm)
}

get_position_probability_matrix <- function(pcm){
    ppm = pcm/colSums(pcm)
    return(ppm)
}

get_background_frequencies_by_gene <- function(data){
    sequences = unique(data$sequence)
    counts = matrix(0, ncol = 1, nrow = 4)
    rownames(counts) = c('A', 'C', 'T', 'G')
    for (base in c('A', 'T', 'C', 'G')){
        counts[base, 1] = 1/(sum(str_count(sequences, base))/sum(nchar(sequences)))
    }
    return(counts)
}

get_log_ppm <- function(data){
    pcm = get_position_count_matrix(data)
    ppm = get_position_probability_matrix(pcm)

    return(log(ppm))
}

get_pwm <- function(data){
    pcm = get_position_count_matrix(data)
    ppm = get_position_probability_matrix(pcm)
    background = get_background_frequencies_by_gene(data)
    background_matrix = matrix(rep(background, 8), 4)
    rownames(background_matrix) = c('A', 'C', 'T', 'G')
    colnames(background_matrix) = paste0('pos', seq(1, 8))
    pwm = log(hadamard.prod(ppm, background_matrix))
    return(pwm)
}

get_p_trim_given_seq <- function(sequence_dataframe){
    denominator = sum(sequence_dataframe$pwm_prob)
    sequence_dataframe$p_trim_given_seq = sequence_dataframe$pwm_prob/denominator
    return(sequence_dataframe)
}

get_trimming_probabilities <- function(sequence, pwm){
    sequence_data = data.table()
    for (trim_length in seq(4, 18)){
        motif = get_motif_from_gene_sequence_trim(sequence, motif_length = 8, trim_length) 
        sequence_data = rbind(sequence_data, data.table(sequence = sequence, trim = trim_length, motif = motif))
    }
    
    sequence_data = separate_motif_column_by_position(sequence_data)
    sequence_data = get_probability_from_pwm(pwm, sequence_data)    
    sequence_data = get_p_trim_given_seq(sequence_data) 
    return(sequence_data)
}

predict_trimming_dist_given_gene <- function(data, model){
    data$predicted_prob = predict(model, type = 'response')
    return(data)
}

get_predicted_dist_file_name <- function(){
    path = file.path(OUTPUT_PATH, 'predicted_trimming_distributions', 'simulated_data')
    dir.create(path)

    complete_path = file.path(path, paste0('simulated_data_predicted_trimming_dist.tsv'))
    return(complete_path)
}


write_predicted_trimming_dist_by_genes <- function(data, genes, model){
    all_genes_predicted_data = predict_trimming_dist_given_gene(data, model)
    file_name = get_predicted_dist_file_name()
    fwrite(all_genes_predicted_data, file_name, sep = '\t')
}

get_predicted_dist_figure_file_name <- function(gene){
    path = file.path(PROJECT_PATH, 'plots', 'predicted_trimming_distributions', 'simulated_data')
    dir.create(path)

    complete_path = paste0(path, '/', gene, '_dist.pdf')
    return(complete_path)
}

get_actual_trimming_dist_by_gene_subject <- function(actual_data, gene_name){
    subset = actual_data[gene == gene_name]
    return(subset)
}

plot_predicted_trimming_dists <- function(actual_data, predicted_data, gene_name){
    actual_dist = get_actual_trimming_dist_by_gene_subject(actual_data, gene_name)
    actual_dist = actual_dist[order(subject, trim)]

    important_cols = c('trim', 'predicted_prob', 'gene')
    predicted_data = predicted_data[gene == gene_name, ..important_cols]
    N = paste(unique(predicted_data$gene_count))
    title = paste0(gene_name, ', N = ', N)
    plot = ggplot() +
        geom_line(data = actual_dist, aes(x = trim, y = actual_p_trim_given_seq, group = subject), size = 2, alpha = 0.7, color = 'grey') +
        geom_line(data = actual_dist, aes(x = trim, y = actual_p_trim_given_seq, group = subject), size = 2, alpha = 0.7) +
        geom_line(data = predicted_data, aes(x = trim, y = predicted_prob), size = 3, alpha = 0.7, color = 'blue') +
        ggtitle(title) +
        xlab('Number of trimmed nucleotides') +
        ylab('Probability') +
        theme_cowplot(font_family = 'Arial') + 
        theme(legend.position = "none", text = element_text(size = 30), axis.text.x=element_text(size = 20), axis.text.y = element_text(size = 20), axis.line = element_blank(),axis.ticks = element_line(color = 'gray60', size = 1.5)) + 
        background_grid(major = 'xy') + 
        panel_border(color = 'gray60', size = 1.5) 
    new_gene_name = str_replace(gene_name, '/', '_')
    plot_name = get_predicted_dist_figure_file_name(new_gene_name) 
    ggsave(plot_name, plot = plot, width = 14, height = 10, units = 'in', dpi = 750, device = cairo_pdf)
    print(paste0('finished plot for ', gene_name))
    return(plot)
}
    
relevel_contrasts <- function(data, ref ){
    for (position in paste0('pos', seq(1, 8))){
        data[[position]] = as.factor(data[[position]])
        data[[position]] = relevel(data[[position]], ref)
    }
    return(data)
}

get_levels <- function(data, ref){
    options(contrasts = rep ("contr.sum", 2)) 
    data = relevel_contrasts(data, ref)
    contrasts = contrasts(data$pos1)
    levels = data.table(base = rownames(contrasts(data$pos1)), number = c(1, 2, 3, NA))
    return(levels)
}

get_coeffiecient_matrix <- function(data, ref){
    data = relevel_contrasts(data, ref)
    model = mclogit(formula, data = data, weights = p_gene * 100, start = rep(0, 24), contrasts = list(pos1 = "contr.sum", pos2 = "contr.sum", pos3 = "contr.sum", pos4 = "contr.sum", pos5 = "contr.sum", pos6 = "contr.sum", pos7 = "contr.sum", pos8 = "contr.sum"))

    levels = get_levels(data, ref)

    together = matrix(0, nrow = 4, ncol = 8)
    colnames(together) = paste0('pos', seq(1, 8))
    rownames(together) = c('A', 'C', 'T', 'G')

    for (coef in seq(1, length(coef(model)))){
        name = names(coef(model)[coef])
        position = substring(name, 1, 4)
        num = substring(name, 5, 5)
        base = levels[number == num]$base
        together[base, position] = coef(model)[coef]
    }

    missing_base = levels[is.na(number)]$base
    together[missing_base,] = -1*colSums(together)

    return(together)
}

