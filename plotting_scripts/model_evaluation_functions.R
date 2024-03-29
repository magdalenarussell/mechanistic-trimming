make_loss_type_names_neat <- function(type) {
    if (type == 'expected_log_loss'){
        nice = 'Expected\nlog loss\n(held-out\ndatasets)'
    } else if (type == 'log_loss'){
        nice = 'Log loss\n(full training\ndataset)'
    } else if (type == 'aic'){
        nice = 'AIC' 
    } else if (type == 'raw_loss'){
        nice = 'Raw count log loss (training dataset)' 
    } else if (type == 'old_loss_cv'){
        nice = 'Raw count conditional log loss\n(held-out dataset)' 
    } else if (type == 'log_loss_j_gene'){
        nice = 'Log loss\n(full J-gene\ndataset)' 
    } else if (type %like% 'v_gene_family_loss'){
        if (type == 'v_gene_family_loss'){
            nice = 'Log loss (\"most-different\"\nheld-out dataset)'
        } else {
            if (type %like% 'full_v_gene_family_loss') {
                nice = str_replace(type, 'full_v_gene_family_loss,', 'Log loss\n(full seq.\nV-gene family\n')
            } else {
                nice = str_replace(type, 'v_gene_family_loss,', 'Log loss\n(V-gene family\n')
            }
            nice = paste0(nice, '\nheld-out)')
        }
    } else {
        nice = type
    }
    return(nice)
}

order_losses <- function(loss_list){
    loss_list_ordered = c(loss_list[loss_list %like% 'training'], 
                          loss_list[loss_list %like% 'Expected'],
                          loss_list[(loss_list %like% 'V-gene family' & !(loss_list %like% 'full seq.')) | (loss_list %like% 'terminal')],
                          loss_list[loss_list %like% 'full seq.'], 
                          loss_list[loss_list %like% 'testing'],
                          loss_list[loss_list %like% 'J-gene'])
    return(unique(loss_list_ordered))
}
 
make_model_names_neat <- function(model_names){
    model_names = str_replace_all(model_names, 'two_side', 'two-side')
    two_side_original = sapply(strsplit(model_names, 'two-side'), `[`, 2)
    two_side_original_subset = unique(two_side_original[!is.na(two_side_original) & !(two_side_original %like% 'base-count') & !(two_side_original %like% 'dinuc-count')])
    two_side_nice = str_replace_all(two_side_original, '_', ' ') 
    two_side_nice_subset = unique(two_side_nice[!is.na(two_side_nice) & !(two_side_nice %like% 'dinuc-count') & !(two_side_nice %like% 'base-count')]) 
    model_names = str_replace_all(model_names, 'dna_shape', 'dna-shape') 
    if (length(two_side_original_subset) > 0) {
        model_names = str_replace_all(model_names, two_side_original_subset, two_side_nice_subset) 
    }
    nice_model_names = str_replace_all(model_names, '_', ' + ')
    return(nice_model_names)
}

make_hairpin_names_neat <- function(motif_types){
    pnucs = sapply(strsplit(motif_types, '_'), `[`, 2) 
    pnucs = str_replace_all(pnucs, '^1', '+1')
    pnucs = str_replace_all(pnucs, '^3', '+3')
    pnucs = paste0(pnucs, ' hairpin opening position')
    pnucs = str_replace_all(pnucs, 'no', 'blunt')
    pnucs= str_replace_all(pnucs, 'NA', '+2')
    return(pnucs)
}

process_model_evaluation_file <- function(eval_data, model_types_neat, left_motif_size_filter = NA, right_motif_size_filter = NA, terminal_melting_5_end_length_filter = NA, lower_trim_bound = LOWER_TRIM_BOUND, upper_trim_bound = UPPER_TRIM_BOUND){
    if (length(model_types_neat) == 1){ 
        if ((model_types_neat %like% 'motif') | (model_types_neat %like% 'shape')){
            sub_model = str_remove(model_types_neat, 'motif_')
            sub_model = str_remove(sub_model, 'dna_shape')
            sub_model = str_remove(sub_model, '-std')
            sub_model = str_replace(sub_model, '__', '_')
            if (substring(sub_model, nchar(sub_model), nchar(sub_model)) == '_'){
                sub_model = substring(sub_model, 1, nchar(sub_model)-1)
            }
            if (substring(sub_model, 1, 1) == '_'){
                sub_model = substring(sub_model, 2, nchar(sub_model))
            }
            eval_data = eval_data[(model_type %in% model_types_neat) | (model_type == sub_model & motif_length_3_end == 0)]
        } else {
            eval_data = eval_data[model_type %in% model_types_neat]
        }
    } else {
        eval_data = eval_data[model_type %in% model_types_neat]
    }

    if (!is.na(left_motif_size_filter)){
        # eval_data = eval_data[(motif_length_5_end == left_motif_size_filter) | (model_type %in% c('distance', 'two_side_terminal_melting', 'distance_two_side_terminal_melting') & motif_length_5_end == 0)]
        eval_data = eval_data[(motif_length_5_end == left_motif_size_filter) | (!((model_type %like% 'motif') | (model_type %like% 'shape')) & motif_length_5_end == 0)]
    }

    if (!is.na(right_motif_size_filter)){
        eval_data = eval_data[(motif_length_3_end == right_motif_size_filter) | (model_type %in% model_types_neat[!((model_types_neat %like% 'motif') | (model_types_neat %like% 'shape'))] & motif_length_3_end ==  0)]
    }

    eval_data = eval_data[terminal_melting_5_end_length %in% terminal_melting_5_end_length_filter]

    eval_data = unique(eval_data)
    eval_data = eval_data[gene_weight_type == GENE_WEIGHT_TYPE]
    eval_data = eval_data[lower_bound %in% lower_trim_bound]
    eval_data = eval_data[upper_bound %in% upper_trim_bound]
    eval_data[loss_type == 'v_gene_family_loss', loss_type := paste0(loss_type, ', cluster ', held_out_clusters)]
    return(eval_data)
}
