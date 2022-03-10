get_base_shape_names <- function(shape, left_base_count, right_base_count){
    if (left_base_count > 0){
        left_bases = paste0(shape, '_5end_pos', seq(left_base_count, 1))    
    } else {
        left_bases = c()
    } 
    if (right_base_count > 0){
        right_bases = paste0(shape, '_3end_pos', seq(1, right_base_count))
    } else {
        right_bases = c()
    }
    positions = c(left_bases, right_bases)
    return(positions)
}

get_bond_shape_names <- function(shape, left_base_count, right_base_count){
    if (left_base_count > 1){
        left_bonds = c(paste0(shape, '_5end_bond', seq(left_base_count-1, 1)))    
    } else {
        left_bonds = c()
    }
    if (left_base_count > 0 & right_base_count > 0) {
        left_bonds = c(left_bonds, paste0(shape, '_trim_bond0'))
    }
    if (right_base_count > 1){
        right_bonds = paste0(shape, '_3end_bond', seq(1, right_base_count-1))
    } else {
        right_bonds = c()
    }
    positions = c(left_bonds, right_bonds)
    return(positions)
}

get_dna_shape_names_positions <- function(shape, left_base_count = LEFT_NUC_MOTIF_COUNT, right_base_count = RIGHT_NUC_MOTIF_COUNT){
    if (shape %in% c('HelT', 'Roll')) {
        positions = get_bond_shape_names(shape, left_base_count, right_base_count)
    } else {
        positions = get_base_shape_names(shape, left_base_count, right_base_count)
    }
    return(positions)
}

subset_dna_shape_data <- function(compiled_dna_shape_data, shapes) {
    positions = c()
    for (shape in shapes){
        shape_positions = get_dna_shape_names_positions(shape)
        positions = c(positions, shape_positions)
    }
    cols = c('window', positions)
    subset = compiled_dna_shape_data[, ..cols]
    return(subset)
}

reformat_dna_shape_data <- function(dna_shape_data, dna_strings, left_window_size, right_window_size){
    window = as.data.table(dna_strings)
    colnames(window) = 'window'
    dna_shape_dt = window
    for (name in names(dna_shape_data)){
        temp = data.table(dna_shape_data[[name]])
        colnames(temp) = get_dna_shape_names_positions(name, left_window_size, right_window_size)
        temp_tog = cbind(temp, window)
        dna_shape_dt = merge(dna_shape_dt, temp_tog, by = 'window')
    }
    
    # subset columns to valid, motif positions
    dna_shape_subset = subset_dna_shape_data(dna_shape_dt, names(dna_shape_data))

    return(dna_shape_subset)
}

get_dna_structure_windows <- function(window_size, left_window_size, right_window_size){
    require(DNAshapeR)    
    require(gtools)
    stopifnot(sum(left_window_size, right_window_size) == window_size)     
    all_combos = permutations(4, window_size, c('A', 'T', 'C', 'G'), repeats.allowed= TRUE) 
    all_combos_strings = apply(all_combos, 1, paste, collapse="")

    temp_dir = paste0(OUTPUT_PATH, '/temp')
    dir.create(temp_dir)
    temp_file_name = paste0(temp_dir, '/dna_structure_windows_size_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '.fa')
    fwrite(as.data.table(all_combos_strings), temp_file_name, col.names=FALSE)

    # predict DNA shape
    dna_shape = getShape(temp_file_name)

    # format DNA shape parameters
    formatted_dna_shape_data = reformat_dna_shape_data(dna_shape, all_combos_strings, left_window_size, right_window_size)

    file_name = paste0(OUTPUT_PATH, '/temp/dna_structure_windows_size_',LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '.tsv')
    fwrite(formatted_dna_shape_data, file_name, sep = '\t')

    # delete excess files
    delfiles = dir(path=temp_dir ,pattern=paste0('dna_structure_windows_size_', LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT, '.fa*'))
    file.remove(file.path(temp_dir, delfiles))
}

process_for_dna_structure <- function(group_motif_data){
    left_window_size = LEFT_NUC_MOTIF_COUNT + 2
    right_window_size = RIGHT_NUC_MOTIF_COUNT + 2
    window_size = left_window_size + right_window_size

    file_name = paste0(OUTPUT_PATH, '/temp/dna_structure_windows_size_',LEFT_NUC_MOTIF_COUNT, '_', RIGHT_NUC_MOTIF_COUNT,'.tsv')

    if (!file.exists(file_name)) {
        print('compiling dna_shape data, first')
        get_dna_structure_windows(window_size, left_window_size, right_window_size)
    }  

    window_shape_data = fread(file_name)

    # merge
    together = merge(group_motif_data, window_shape_data, by.x = 'motif', by.y = 'window')

    return(together)
}
