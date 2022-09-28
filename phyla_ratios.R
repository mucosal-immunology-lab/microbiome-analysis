# Function to calculate the ratio of different phyla
phyla_ratios <- function(phyloseq_object = NULL, test_variable = NULL, average_reads_threshold = 1000,
                         plot_output_folder = NULL, plot_file_prefix = NULL) {
  # Load required packages
  library(pacman)
  pkgs <- c('ggplot2', 'IRanges', 'base', 'BiocGenerics', 'Matrix', 'S4Vectors', 'biomformat', 'rstatix', 
            'ggpubr', 'plotly', 'dplyr', 'stats', 'phyloseq', 'ggbeeswarm', 'ggpmisc', 'ggsci', 'grDevices', 
            'here', 'purrr', 'stringr', 'tibble', 'utils')
  pacman::p_load(char = pkgs)
  
  # Function to rotate a data.frame and maintain names
  rotate_df <- function(data){
    names_row <- rownames(data)
    names_col <- colnames(data)
    
    data_rotated <- data.frame(t(data))
    rownames(data_rotated) <- names_col
    colnames(data_rotated) <- names_row
    
    data_rotated
  }
  
  # Extract tax table data
  tax_table <- data.frame(tax_table(phyloseq_object))
  
  # Extract test variable from the sample data
  if (!is.null(test_variable)) {
    test_var <- sample_data(phyloseq_object)[[test_variable]]
  } else {
    stop('Please provide some grouping variable to the "test_variable" argument.\n
         This must be the name of the column from the sample_data element of the phyloseq object,
         and can be of class factor, character, or numeric.')
  }
  
  # Detect the phylum column name and agglomerate to that level
  phylum_colname <- colnames(tax_table)[str_detect(colnames(tax_table), 'hylum')]
  if (length(phylum_colname) == 0) {
    stop('Could not detect a tax_table column name that matched either "Phylum" or "phylum".')
  }
  physeq_phylum <- phyloseq::tax_glom(phyloseq_object, taxrank = phylum_colname)
  
  # Reextract the tax table and otu table data
  otu_table <- data.frame(otu_table(physeq_phylum))
  rownames(otu_table) <- tax_table(physeq_phylum)[, phylum_colname]
  otu_table <- otu_table[!str_detect(rownames(otu_table), 'nknown'), ] # remove Unknowns
  
  # Add a rowSums column, and filter by the average_reads_threshold
  reads_threshold <- average_reads_threshold * ncol(otu_table) # get minimum rowSums by multiplying columns by average_reads_threshold
  otu_table <- otu_table %>%
    mutate(row_sum = rowSums(otu_table)) %>%
    filter(row_sum >= reads_threshold) %>%
    dplyr::select(-row_sum) %>%
    rotate_df()
  
  # Prepare all combinations of the phyla
  combinations <- flatten(lapply(seq_along(colnames(otu_table)), 
                                 function(x) combn(colnames(otu_table), x, FUN = list)))
  combinations_len2 <- combinations[sapply(combinations, length) == 2]
  
  # Prepare a list of subset data.frames containing just the two columns (plus reverse directions)
  phyla_comparisons <- list()
  for (i in seq_along(combinations_len2)) {
    pair <- combinations_len2[[i]]
    pair_name <- paste(pair, collapse = '/')
    pair_rev <- rev(pair)
    pair_rev_name <- paste(pair_rev, collapse = '/')
    
    pair_otu_table <- dplyr::select(otu_table, pair) %>%
      mutate(ratio = .data[[pair[1]]] / .data[[pair[2]]],
             test_var = test_var)
    
    pair_rev_otu_table <- dplyr::select(otu_table, pair_rev) %>%
      mutate(ratio = .data[[pair_rev[1]]] / .data[[pair_rev[2]]],
             test_var = test_var)
    
    phyla_comparisons[[pair_name]] <- pair_otu_table
    phyla_comparisons[[pair_rev_name]] <- pair_rev_otu_table
  }
  
  # Prepare plots for each of the phyla comparisons (including reverse directions)
  phyla_comparison_plots <- list()
  phyla_comparison_stats <- list()
  for (i in seq_along(names(phyla_comparisons))) {
    comp <- names(phyla_comparisons)[i]
    comp_df <- phyla_comparisons[[comp]]
    
    # Run stats
    comp_wilcox <- wilcox_test(comp_df, ratio ~ test_var)
    phyla_comparison_stats[[comp]] <- comp_wilcox
    
    # Set utility variables
    p_ydist <- diff(range(comp_df$ratio))
    title_lab <- paste0(comp, ' Ratio')
    
    if (class(comp_df$test_var) %in% c('character', 'factor')) {
      # Plot
      p <- ggplot(comp_df, aes(x = test_var, y = ratio)) +
        geom_boxplot(aes(fill = test_var)) +
        geom_beeswarm(binaxis = 'y', stackdir = 'center', binwidth = p_ydist / 30, cex = 3) +
        scale_fill_jama(name = test_variable, alpha = 0.6) +
        guides(fill = 'none') +
        stat_compare_means(size = 2) +
        labs(title = title_lab,
             x = test_variable,
             y = 'Ratio') +
        theme(text = element_text(size = 8))
      
      phyla_comparison_plots[[comp]] <- p
    }
    
    if (class(comp_df$test_var) == 'numeric') {
      # Check direction of association
      temp_lm <- lm(comp_df$ratio ~ comp_df$test_var)
      direction = ifelse(temp_lm$coefficients[['comp_df$test_var']] > 0, 'blue', 'red')
      
      # Plot
      p <- ggplot(comp_df, aes(x = test_var, y = ratio)) +
        geom_smooth(method = 'lm', formula = y ~ x, se = FALSE, color = direction) +
        stat_poly_eq(method = 'lm', formula = y ~ x, size = 2) +
        geom_point() +
        guides(color = 'none') +
        labs(title = title_lab,
             x = test_variable,
             y = 'Ratio') +
        theme(text = element_text(size = 8))
      
      phyla_comparison_plots[[comp]] <- p
    }
  }
  
  # Combine the stats output
  stats_combined <- do.call(rbind, phyla_comparison_stats) %>%
    rownames_to_column(var = 'comparison') %>%
    mutate(comparison = gsub('\\.\\d*', '', comparison))
  stats_combined_signif <- stats_combined %>%
    filter(p < 0.05)
  
  # Save plots
  if (!is.null(plot_output_folder)) {
    if (!is.null(plot_file_prefix)) {
      plot_fp <- here::here(plot_output_folder, paste(plot_file_prefix, 'phyla_ratios.pdf', sep = '_'))
    } else {
      plot_fp <- here::here(plot_output_folder, 'phyla_ratios.pdf')
    }
    
    if (length(phyla_comparison_plots) > 0) {
      plot_pages <- ggarrange(plotlist = phyla_comparison_plots, nrow = 4, ncol = 3)
      
      pdf(file = plot_fp, width = 8.5, height = 11, bg = 'white')
      if (class(plot_pages)[1] != 'list') {
        print(plot_pages)
      } else {
        for (i in 1:length(plot_pages)) {
          print(plot_pages[[i]])
        }
      }
      dev.off()
    }
  }
  
  # Return object
  return(list(phyla_ratios = phyla_comparisons,
              phyla_ratio_plots = phyla_comparison_plots,
              stats_all = stats_combined,
              stats_signif = stats_combined_signif))
}