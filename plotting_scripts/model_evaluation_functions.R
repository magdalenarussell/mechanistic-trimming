filter_model_types <- function(remove_types_with_string = NA){
    model_type_files = list.files(path = 'scripts/model_formula_functions/')
    model_types = str_sub(model_type_files[model_type_files != '_ignore' & model_type_files != "model_formula_specific_functions"], end = -3)
    model_types_neat = model_types
    if (!(is.na(remove_types_with_string))){
        for (string in remove_types_with_string){
            if (string == 'terminal_melting'){
                model_types_neat = model_types_neat[!grepl(string, model_types_neat) | grepl('two_side_terminal_melting', model_types_neat)]
            } else {
                model_types_neat = model_types_neat[!grepl(string, model_types_neat)]
            }
        }
    }
    model_types_neat = model_types_neat[order(model_types_neat)]
    return(model_types_neat)
}

process_model_evaluation_file <- function(eval_data, model_types_neat, left_motif_size_filter = NA, right_motif_size_filter = NA, terminal_melting_5_end_length_filter = NA){
    if (length(model_types_neat) == 1){ 
        if (model_types_neat %like% 'motif'){
            sub_model = str_split(model_types_neat, 'motif_')[[1]][2]
            eval_data = eval_data[(model_type %in% model_types_neat) | (model_type == sub_model & motif_length_3_end == 0)]
        } else {
            eval_data = eval_data[model_type %in% model_types_neat]
        }
    } else {
        eval_data = eval_data[model_type %in% model_types_neat]
    }

    if (!is.na(left_motif_size_filter)){
        eval_data = eval_data[(motif_length_5_end == left_motif_size_filter) | (model_type %in% c('distance', 'two_side_terminal_melting', 'distance_two_side_terminal_melting') & motif_length_5_end == 0)]
    }

    if (!is.na(right_motif_size_filter)){
        eval_data = eval_data[(motif_length_3_end == right_motif_size_filter) | (model_type %in% model_types_neat[!(model_types_neat %like% 'motif')] & motif_length_3_end ==  0)]
    }

    eval_data = eval_data[terminal_melting_5_end_length %in% terminal_melting_5_end_length_filter]

    eval_data = unique(eval_data)
    eval_data = eval_data[gene_weight_type == GENE_WEIGHT_TYPE]
    eval_data = eval_data[lower_bound == LOWER_TRIM_BOUND]
    eval_data = eval_data[upper_bound == UPPER_TRIM_BOUND]
    return(eval_data)
}

get_term_count <- function(processed_eval_data, right_motif_count, left_motif_count){
    total_motif_size = right_motif_count + left_motif_count
    number_of_distance_terms = UPPER_TRIM_BOUND - LOWER_TRIM_BOUND
    processed_eval_data[, terms := 0]
    processed_eval_data[model_type %like% 'motif', terms := terms + total_motif_size*4]
    processed_eval_data[model_type %like% 'distance', terms := terms + number_of_distance_terms]
    processed_eval_data[model_type %like% 'terminal_melting', terms := terms + 1]
    processed_eval_data[model_type %like% 'terminal_gc', terms := terms + 1]
    processed_eval_data[model_type %like% 'two_side_terminal_melting', terms := terms + 1]
    processed_eval_data[model_type %like% 'two_side_base_count', terms := terms + 8]
    return(processed_eval_data)
}

get_terminal_melting_calculation_type <- function(processed_eval_data){
    processed_eval_data[model_type %like% 'terminal_melting_', melting_type := sapply(model_type, function(x) tail(str_split(x, '_')[[1]], 1))]
    processed_eval_data[model_type %like% 'terminal_melting' & is.na(melting_type), melting_type := 'simple']
    processed_eval_data[is.na(melting_type), melting_type := 'NA']

    processed_eval_data[model_type %like% 'terminal_melting_NN', model_type := substring(model_type, 1, nchar(model_type)- 3)]
    processed_eval_data[model_type %like% 'terminal_melting_combo', model_type := substring(model_type, 1, nchar(model_type)- 6)]
    return(processed_eval_data)
}

get_pnuc_count <- function(processed_eval_data){
    processed_eval_data[motif_type != 'unbounded' & motif_type != 'unbounded_no_pnuc', pnuc_count := as.numeric(sapply(motif_type, function(x) str_split(x, '_')[[1]][2]))]
    processed_eval_data[motif_type != 'unbounded' & motif_type != 'unbounded_no_pnuc' & pnuc_count > 0, nick_position := paste0('+', pnuc_count)]
    processed_eval_data[motif_type != 'unbounded' & motif_type != 'unbounded_no_pnuc' & pnuc_count < 0, nick_position := paste0('-', pnuc_count)]
    processed_eval_data[motif_type == 'unbounded_no_pnuc', c("pnuc_count", "nick_position"):= list(0, '0')]
    processed_eval_data[motif_type == 'unbounded', c("pnuc_count", "nick_position") := list(2, '+2')]
    return(processed_eval_data)
}
