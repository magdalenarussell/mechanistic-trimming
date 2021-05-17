get_PWM <- function(conditioning){
    file = get_fitted_PWM_path(conditioning)
    PWM = as.matrix(fread(file), rownames = 1)
    return(PWM)
}


create_plot_title <- function(conditioning){
    if (REGRESSION_TYPE == 'all_subject'){
    title = paste0(PFM_TYPE, ' ', MOTIF_TYPE, ' conditioned on ', conditioning, ' for ', REGRESSION_TYPE)
    } 
    #TODO add per subject stuff
    return(title)
}

create_file_name <- function(conditioning){
    name = paste0(paste0('PWM_', REGRESSION_TYPE, '_', MOTIF_TYPE, '_', PFM_TYPE, '_', conditioning, '_conditioning.pdf'))
    path = file.path('plots', name)
    return(path)
}

remap_position_labels <- function(PWM){
    seq = c(seq(LEFT_NUC_MOTIF_COUNT, 1), seq(-1, -RIGHT_NUC_MOTIF_COUNT))
    colnames(PWM) = seq
    return(PWM)
}

transform_PWM_to_dataframe <- function(PWM){
    PWM = remap_position_labels(PWM)
    df = PWM %>%
        as.data.frame() %>%
        rownames_to_column('Base') %>%
        pivot_longer(-c(Base), names_to = "Position", values_to = "weight") 
    return(df)
}


plot_single_PWM_heatmap <- function(conditioning, subtitle = NULL){
    PWM = get_PWM(conditioning)
    plot_title = create_plot_title(conditioning)
    filename = create_file_name(conditioning)
    
    weighted_PWM = PWM/log(10)
    PWM_df = transform_PWM_to_dataframe(weighted_PWM)
    
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

#TODO update the following!!
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

