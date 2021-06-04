library(cowplot)
library(ggplot2)
library(data.table)
library(tidyverse)
library(mclogit)
library(matrixcalc)
source('config/config.R')

source('create_simulation_data/simulation_functions.R')

# First, generate random sequences and calculate empirical PWM
set.seed(60)
trims = seq(4, 18)

all_subjects = data.table()
for (gene in paste0('gene_', seq(1, 10))){
    sequence = rep(get_random_gene_sequence(8, 4, 18), 16)

    motifs = c()
    for (index in seq(length(trims))){
        motif = get_motif_from_gene_sequence_trim(sequence[index], 8, trims[index])    
        motifs = c(motifs, motif)
    }

    count = sample(1000:2000, 16, replace=T)
    
    together = data.table(gene = rep(gene, 16), trim_length = trims, sequence = sequence, motif = motifs, count = count)
    all_subjects = rbind(all_subjects, together)
}

final = separate_motif_column_by_position(all_subjects)
final = create_position_base_variables(final)

PWM = get_pwm(final)


# Second, with this PWM, generate data (i.e. draw a gene sequence and identify probability of trimming at each trim_length using PWM)
set.seed(111)
all_genes = data.table()
for (gene in paste0('gene_', seq(1, 60))){
    sequence = get_random_gene_sequence(8, 4, 18)
    trimming_probs = get_trimming_probabilities(sequence, PWM)
    for (subject in paste0('subject', seq(1, 100))){
        # Now, draw a gene count (number of trimming events)
        trimming_event_count = sample(1000:2000, 1, replace=T)
    
        # Now with these probabilities, simulate trimming events
        trimming_events = as.data.table(table(sample(trimming_probs$trim, trimming_event_count, prob = trimming_probs$p_trim_given_seq, replace = TRUE)))
        colnames(trimming_events) = c('trim', 'trim_count')
        trimming_events$trim = as.integer(trimming_events$trim)
        
        together = merge(trimming_probs, trimming_events, fill = TRUE)
        together$gene = gene
        together$subject = subject
        together$gene_count = trimming_event_count
        together$actual_p_trim_given_seq = together$trim_count/together$gene_count

        all_genes = rbind(all_genes, together)
    }
}

all_genes[, p_gene := gene_count/sum(gene_count), by = .(trim, subject)]
all_genes = create_position_base_variables(all_genes)

eq_string = paste0('pos', seq(1, 8))
eq = paste(eq_string, collapse = ' + ')
formula = formula(paste0('cbind(trim_count, interaction(gene, subject)) ~ ', eq))

model = mclogit(formula, data = all_genes, weights = p_gene * 100, start = rep(0, 24), contrasts = list(pos1 = "contr.sum", pos2 = "contr.sum", pos3 = "contr.sum", pos4 = "contr.sum", pos5 = "contr.sum", pos6 = "contr.sum", pos7 = "contr.sum", pos8 = "contr.sum"))

predicted_data = predict_trimming_dist_given_gene(all_genes, model)
for (gene_name in unique(all_genes$gene)){
    plot_predicted_trimming_dists(all_genes, predicted_data, gene_name)
}
