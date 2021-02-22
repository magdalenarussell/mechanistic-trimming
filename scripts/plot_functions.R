get_PFM <- function(group_name, weighting = NULL){
    output_path = file.path('_ignore', 'pfms', 'complete') 
    file = file.path(output_path, paste0(group_name, '_', PFM_TYPE, '_', MOTIF_TYPE, '_', weighting, '.tsv'))
    PFM = as.matrix(fread(file), rownames = 1)
    return(PFM)
}


create_plot_title <- function(group_name, weighting = NULL){
    if (group_name == 'ALL'){
        title = paste0(PFM_TYPE, ' ', MOTIF_TYPE, ' for ', group_name, ' subjects')
    } else{
        snp = strsplit(group_name, '_')[[1]][1]
        title = paste0(PFM_TYPE, ' ', MOTIF_TYPE, ' by genotype for SNP ', snp)
    }
    title = ifelse(is.null(weighting), title, paste0(title, ' with gene/nt weighting'))
    return(title)
}

create_file_name <- function(group_name, weighting = NULL, plot_type){
    output_path = file.path('plots') 
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
    filename = file.path(output_path, paste0(group_name, '_', PFM_TYPE, '_', MOTIF_TYPE, '_', weighting, '_', plot_type, '.pdf'))
    return(filename)
}

transform_PWM_to_dataframe <- function(PWM){
    df = PWM %>%
        as.data.frame() %>%
        rownames_to_column('Base') %>%
        pivot_longer(-c(Base), names_to = "Position", values_to = "weight") 
    return(df)
}

plot_single_PWM_heatmap <- function(group_name, weighting = NULL, subtitle = NULL){
    PFM = get_PFM(group_name, weighting)
    plot_title = create_plot_title(group_name, weighting)
    filename = create_file_name(group_name, weighting, 'PWM_heatmap_with_background_correction')
    
    PWM = calculate_PWM_from_PFM(PFM, weighting)/log(10)
    PWM_df = transform_PWM_to_dataframe(PWM)
    
    PWM_df$Base = factor(PWM_df$Base, levels = c('T', 'G', 'C', 'A'))
    PWM_df$Position = factor(PWM_df$Position, levels = sort(unique(as.numeric(PWM_df$Position))))

    PWM_df %>%
        ggplot(aes(x=Position, y=Base, fill=weight)) + 
        geom_tile() +
        theme_classic() +
        labs(title = plot_title, subtitle = subtitle) +
        theme(text = element_text(size = 20)) +
        geom_vline(xintercept = LEFT_NUC_MOTIF_COUNT + 0.5, size = 3, color = 'black') +
        # ggtitle(plot_title) +
        guides(fill = guide_colourbar(barheight = 10)) +
        scale_fill_viridis_c(name = 'log10(probability of deletion)')
    ggsave(filename, plot = last_plot(), width = 10, height = 5, units = 'in', dpi = 750, device = 'pdf')
}

plot_PWM_heatmap_by_genotype <- function(snpID, weighting = NULL, subtitle = NULL){
    together = data.frame()
    for (genotype in c(0,1,2)){
        assign(paste0('PFM', genotype), get_PFM(group_name = get_group_name(snpID, genotype) , weighting))
        assign(paste0('PWM', genotype), calculate_PWM_from_PFM(get(paste0('PFM', genotype)), weighting)/log(10))
        PWM_df = transform_PWM_to_dataframe(get(paste0('PWM', genotype)))
        PWM_df$genotype = genotype
        together = rbind(together, PWM_df)
    }

    plot_title = create_plot_title(group_name = paste(snpID), weighting)
    
    filename = create_file_name(group_name = paste(snpID), weighting, 'PWM_heatmap_by_GENOTYPE_with_background_correction')

    together$Base = factor(together$Base, levels = c('T', 'G', 'C', 'A'))
    together$Position = factor(together$Position, levels = sort(unique(as.numeric(together$Position))))

    together %>%
        arrange(Base) %>%
        ggplot(aes(x=Position, y=Base, fill=weight)) + 
        facet_wrap(~genotype, ncol = 1) +
        geom_tile() +
        theme_classic() +
        labs(title = plot_title, subtitle = subtitle) +
        theme(text = element_text(size = 20)) +
        geom_vline(xintercept = LEFT_NUC_MOTIF_COUNT + 0.5, size = 3, color = 'black') +
        # ggtitle(plot_title) +
        guides(fill = guide_colourbar(barheight = 10)) +
        scale_fill_viridis_c(name = 'log10(probability of deletion)')
    ggsave(filename, plot = last_plot(), width = 10, height = 15, units = 'in', dpi = 750, device = 'pdf')
}

generate_subtitle_from_GWAS_pvalue_slope <- function(snpID, GWAS_results){
    subset = GWAS_results[snp == snpID]
    subtitle = paste0()
    for (row in seq(1, nrow(subset))){
        pvalue = paste0('p = ', signif(subset[row,]$pvalue, digits = 3))
        slope = paste0('slope = ', signif(subset[row,]$slope, digits = 3))
        productivity = ifelse(subset[row,]$productive == 'FALSE', 'non-productive', 'productive')
        phenotype = paste0(' for ', productivity, ' ', subset[row,]$phenotype)
        subtitle = paste0(subtitle, '\n', pvalue, ', ', slope, phenotype)
    }
    return(subtitle)
}

