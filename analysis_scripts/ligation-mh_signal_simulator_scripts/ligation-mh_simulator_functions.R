get_igor_gene_usage_params <- function(){
    # get igor params
    j_choice_path= paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_jchoice_params.tsv')
    v_choice_path= paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_vchoice_params.tsv')

    jg = fread(j_choice_path)
    vg = fread(v_choice_path)
    return(list(v_choice = vg, j_choice = jg))
}

get_trimming_probs <- function(type){
    stopifnot(type %in% c('igor', 'motif_two-side-base-count-beyond'))
    if (type == 'igor'){
        jtrim_path = paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_j_trim_params.tsv')
        vtrim_path = paste0(MOD_OUTPUT_PATH, '/meta_data/igor_alpha_v_trim_params.tsv')

        jt = fread(jtrim_path)
        vt = fread(vtrim_path)
    } else {
        # get predictions (i.e. probabilities) from no-MH model
        v_probs_path = get_model_predictions_file_path('False', model_type='motif_two-side-base-count-beyond')
        j_probs_path = get_model_predictions_file_path('False', model_type='motif_two-side-base-count-beyond')
        v_probs_path = str_replace(v_probs_path, ANNOTATION_TYPE, 'igor_sim_alpha')
        j_probs_path = str_replace(j_probs_path, ANNOTATION_TYPE, 'igor_sim_alpha')
        v_probs_path = str_replace(v_probs_path, PARAM_GROUP, 'nonproductive_v_trim')
        j_probs_path = str_replace(j_probs_path, PARAM_GROUP, 'nonproductive_j_trim')

        if (!file.exists(v_probs_path)){
            print(paste0('need to produce predictions file: ', v_probs_path, ' before trimming probabilities can be obtained'))
        }

        if (!file.exists(j_probs_path)){
            print(paste0('need to produce predictions file: ', j_probs_path, ' before trimming probabilities can be obtained'))
        }

        vt = fread(v_probs_path)
        vcols = c('v_gene', 'v_trim', 'predicted_prob')
        vt = unique(vt[, ..vcols])
        setnames(vt, 'predicted_prob', 'v_trim_prob')

        jt = fread(j_probs_path)
        jcols = c('j_gene', 'j_trim', 'predicted_prob')
        jt = unique(jt[, ..jcols])
        setnames(jt, 'predicted_prob', 'j_trim_prob')
    }
    vt = vt[v_trim >= LOWER_TRIM_BOUND & v_trim <= UPPER_TRIM_BOUND]
    jt = jt[j_trim >= LOWER_TRIM_BOUND & j_trim <= UPPER_TRIM_BOUND]
    return(list(v_trim = vt, j_trim = jt))
}

get_ligation_probabilities <- function(lig_param, possible_ligs){
    stopifnot(lig_param >= 0)
    intercept = 1
    lig_probs = possible_ligs*lig_param + intercept
    names(lig_probs) = possible_ligs    
    lig_probs['no_ligation'] = 0.1
    prob_sum = sum(lig_probs)
    lig_probs = lig_probs/prob_sum
    return(lig_probs)
}

get_all_possible_configs <- function(gene_probs, trim_probs){
    # get all gene and trimming combos
    jcols = c('j_gene', 'j_gene_sequence')
    vcols = c('v_gene', 'v_gene_sequence')
    
    jgene_seq = gene_probs$j_choice[, ..jcols]
    vgene_seq = gene_probs$v_choice[, ..vcols]

    jgene_seq_list = split(jgene_seq, seq(nrow(jgene_seq)))
    vgene_seq_list = split(vgene_seq, seq(nrow(vgene_seq)))

    trims = unique(trim_probs$v_trim$v_trim)

    combos = CJ(vgene_seq_list, jgene_seq_list, v_trim = trims, j_trim = trims, sorted = FALSE)

    combos[, c('v_gene', 'v_gene_sequence') := .(unlist(sapply(vgene_seq_list, function(dt) dt$v_gene)),
                                                 unlist(sapply(vgene_seq_list, function(dt) dt$v_gene_sequence)))]
    combos[, c('j_gene', 'j_gene_sequence') := .(unlist(sapply(jgene_seq_list, function(dt) dt$j_gene)),
                                                 unlist(sapply(jgene_seq_list, function(dt) dt$j_gene_sequence)))]

    cols = colnames(combos)[!(colnames(combos) %like% 'list')]
    combos = combos[, ..cols]

    lig_mat = matrix(0, nrow = nrow(combos), ncol = length(seq(1, 15)))

    for (overlap in seq(ncol(lig_mat))){
        lig = get_possible_ligation_mh_fixed_trim(combos, overlap_count = overlap)
        lig_mat[, overlap] = lig
    }

    # get unique vals 
    unique_mh_list = apply(lig_mat, 1, unique)

    # combine
    possible = data.table(combos, possible_ligation_mh = unique_mh_list)

    # Expand the rows
    expanded_possible = as.data.table(possible[, unnest(.SD, cols = c("possible_ligation_mh"))])

    return(expanded_possible)
}

get_possible_zero_mh_sites <- function(){
    # This function will only return frame data for annotations that are observable, after re-annotation to accommodate cases of maximal ligation-mh
    # Read only necessary columns and apply filter
    frames = fread('https://raw.githubusercontent.com/phbradley/conga/master/conga/tcrdist/db/combo_xcr.tsv', select = c('id', 'region', 'nucseq', 'organism', 'chain', 'frame'))[organism == 'human' & chain == substring(CHAIN_TYPE, 3, 3)]

    # Combine operations to reduce redundancy
    v = frames[region == 'V', .(v_gene = id, v_frame = frame, v_seq = nucseq)]
    j = frames[region == 'J', .(j_gene = id, j_frame = frame, j_seq = nucseq)]

    # Merge operations
    v$dummy = 1
    j$dummy = 1
    gene_pairs = merge(v, j, by = 'dummy', allow.cartesian = TRUE)

    # Get all trimming sites
    trims = data.table(expand.grid(v_trim = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND + 10), 
                                   j_trim = seq(LOWER_TRIM_BOUND, UPPER_TRIM_BOUND + 10)), 
                                   ligation_mh = 0)
    trims$dummy = 1

    # Merge genes and trims
    all = merge(gene_pairs, trims, by = 'dummy', allow.cartesian = TRUE)[, -c('dummy')]
    # Get sequence lengths and subset data to necessary columns
    all[, c('v_seq_len', 'j_seq_len') := .(nchar(v_seq), nchar(j_seq))]
    cols = c('v_gene', 'j_gene', 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'v_trim', 'j_trim', 'ligation_mh')
    all = unique(all[, ..cols])

    # Get oriented full sequences and group genes by common features, also subset columns again
    cols2 = c('v_gene', 'j_gene', 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'v_trim', 'j_trim', 'ligation_mh', 'v_gene_sequence', 'j_gene_sequence')
    adjusted_grouped = get_oriented_full_sequences(all)
    adjusted_grouped = unique(adjusted_grouped[, ..cols2])
    
    # Filter by trimming length
    adjusted_grouped = adjusted_grouped[v_trim <= UPPER_TRIM_BOUND & j_trim <= UPPER_TRIM_BOUND]

    # Get stop positions
    adjusted_grouped = get_stop_codon_positions(adjusted_grouped)

    # Frame calculations
    ## j_frame is subtracted because the overall sequence length should have a count of j_frame excess nucleotides 
    adjusted_grouped[, overall_frame := (v_seq_len + j_seq_len - (v_trim + j_trim) - ligation_mh - j_frame) %% 3]
    adjusted_grouped[, frame_type := fifelse(overall_frame == 0, 'In', 'Out')]

    # Look for stop codons based on frame
    adjusted_grouped = get_stop_positions_with_frame(adjusted_grouped)

    # Subset data by columns
    cols3 = c('v_gene', 'j_gene', 'v_frame', 'j_frame', 'v_seq_len', 'j_seq_len', 'v_trim', 'j_trim', 'ligation_mh', 'overall_frame', 'frame_type', 'frame_stop')
    final = adjusted_grouped[, ..cols3]
    possible_sites = final[frame_type == 'Out' | frame_stop == TRUE] 

    # Define columns for filtering possible sites
    cols = c('v_gene', 'j_gene', 'frame_type', 'overall_frame', 'frame_stop', 'v_trim', 'j_trim', 'ligation_mh')
    possible_sites_subset = unique(possible_sites[, ..cols])
    return(possible_sites_subset)
}
