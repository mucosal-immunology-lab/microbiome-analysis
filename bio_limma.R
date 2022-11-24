############################################################################################
# Copyright (c) 2022 - Mucosal Immunology Lab, Monash University, Melbourne, Australia     #
# Author: Matthew Macowan                                                                  #
# This script is provided under the MIT licence (see LICENSE.txt for details)              #
############################################################################################

bio_limma <- function(input_data, metadata_var = NULL, metadata_condition = NULL, model_matrix = NULL,
                      model_formula_as_string = NULL, use_contrast_matrix = TRUE, coefficients = NULL,
                      factor_reorder_list = NULL, continuous_modifier_list = NULL,
                      contrast_matrix = NULL, adjust_method = 'BH', rownames = NULL, 
                      tax_id_col = NULL, adj_pval_threshold = 0.05, logFC_threshold = 1, 
                      legend_metadata_string = NULL, volc_plot_title = NULL, volc_plot_subtitle = NULL,
                      volc_plot_xlab = NULL, volc_plot_ylab = NULL, remove_low_variance_taxa = FALSE,
                      plot_output_folder = NULL, plot_file_prefix = NULL, redo_boxplot_stats = FALSE) {
  # Load required packages
  pkgs <- c('BiocGenerics', 'base', 'ggtree', 'ggplot2', 'IRanges', 'Matrix', 'S4Vectors', 'biomformat', 'plotly', 
            'dplyr', 'rstatix', 'ggpubr', 'stats', 'phyloseq', 'SummarizedExperiment', 'ggpmisc', 'ggrepel', 'ggsci', 
            'grDevices', 'here', 'limma', 'stringr')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Function to rotate a data.frame and maintain names
  rotate_df <- function(data){
    names_row <- rownames(data)
    names_col <- colnames(data)
    
    data_rotated <- data.frame(t(data))
    rownames(data_rotated) <- names_col
    colnames(data_rotated) <- names_row
    
    data_rotated
  }
  
  # Perform sanity checks
  if (is.null(model_formula_as_string) && is.null(metadata_var)) {
    stop('Please provide at least one of the following: a) preferably a model formula as a string 
         (arg = model_formula_as_string), or b) the name of a single sample_data column (arg = metadata_var).')
  }
  if (!is.null(model_formula_as_string) && !is.null(metadata_var)) {
    metadata_var <- NULL
    message('Please note that because you provided inputs for both model_formula_as_string and metadata_var, 
            only the model_formula_as_string input will be used.')
  }
  if (!is.null(metadata_var) && length(metadata_var) > 1) {
    stop('When using the metadata_var argument, you can only select a single variable from the phyloseq object\'s
         sample_data. If you want to include more, you can provide an input to model_formula_as_string instead.')
  }
  if (!is.null(model_formula_as_string) && is.null(coefficients)) {
    stop('As you have selected to use coefficients instead of a contrast matrix, please assign indices to the 
           coefficients parameter to select coefficients you want to retain for Bayes statistics.')
  }
  if (!dir.exists(plot_output_folder)) {dir.create(plot_output_folder)}
  # Set 'use_contrast_matrix' to FALSE if coefficients are provided
  if (!is.null(coefficients)) {
    use_contrast_matrix <- FALSE
  }
  
  # Define a function to ensure that character 'NA' values are switched with real 'NA' values
  ensure_NA <- function(vector) {
    replace(vector, vector == 'NA', NA)
  }
  
  # Define a function to determine whether any 'NA' values exist on a given row
  any_NA <- function(data) {
    keep_vector <- c()
    for (i in 1:nrow(data)) {
      if (sum(is.na(data[i, ])) > 0) {
        keep_vector <- c(keep_vector, TRUE)
      } else {
        keep_vector <- c(keep_vector, FALSE)
      }
    }
    return(keep_vector)
  }
  
  # Prepare input data depending on input data type
  if (class(input_data) == 'phyloseq') {
    if (nrow(otu_table(input_data)) != nrow(tax_table(input_data))) { # Phyloseq-specific check
      stop('Please ensure that the OTU table and tax table of your phyloseq object have the same number of rows. 
         This should be the case if you have generated your phyloseq object correctly, and ensured that 
         taxa_are_rows = TRUE. I.e. "OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)".')
    }
    
    # Filter phyloseq object by the conditional statement if present
    if(!is.null(metadata_condition)) {
      input_data <- prune_samples(metadata_condition, input_data)
    }
    
    # Check the the tax_id_col is equal to the most deepest classification where not all value are NA where provided
    input_tax_table <- data.frame(tax_table(input_data))
    input_tax_table_not_na <- names(which(colSums(data.frame(is.na(input_tax_table))) < nrow(input_tax_table)))
    max_input_tax_level <- input_tax_table_not_na[length(input_tax_table_not_na)]
    
    if (!is.null(tax_id_col)) {
      if (tax_id_col != max_input_tax_level) {
        stop('The tax ID column you have provided does not match the name of the tax_table column
           with the deepest taxonomic classification for which non-NA values exist.
           Either let the function choose the tax_id_col for you, or if you want to test at a
           different taxonomic level, please use the phyloseq::tax_glom function to agglomerate
           your data.')
      }
    } else if (is.null(tax_id_col)) {
      tax_id_col <- max_input_tax_level
    }
    
    # Prune unknown taxa
    input_data <- prune_taxa(!str_detect(tax_table(input_data)[, tax_id_col], 'Unknown'), input_data)
    
    # Remove low variance taxa where requested
    if (isTRUE(remove_low_variance_taxa)) {
      otu_variance <- apply(data.frame(otu_table(input_data)), 1, var)
      non_zero_variance <- otu_variance > 0
      input_data <- prune_taxa(non_zero_variance, input_data)
    }
    
    # Prepare limma data.frame
    bio_limma_data <- data.frame(otu_table(input_data))
    
    # Set rownames to tax_id_col unless provided another vector
    if (is.null(rownames)) {
      rownames(bio_limma_data) <- make.unique(tax_table(input_data)[,tax_id_col])
    } else {
      rownames(bio_limma_data) <- rownames
    }
    
    # Ensure that NA values are actually NA, and not character type 'NA'
    if (!is.null(metadata_var)) {
      if (!is.null(metadata_condition)) {
        test_df <- data.frame(sample_data(input_data)[, metadata_var])[metadata_condition,]
      } else {
        test_df <- data.frame(sample_data(input_data)[, metadata_var])
      }
    } else {
      if (!is.null(metadata_condition)) {
        test_df <- data.frame(sample_data(input_data))[metadata_condition,]
      } else {
        test_df <- data.frame(sample_data(input_data))
      }
    }
  }
  if (class(input_data) == 'SummarizedExperiment') {
    # Filter assay data by the conditional statement if present
    if(!is.null(metadata_condition)) {
      input_data <- input_data[, metadata_condition]
    } else {
      input_data <- input_data
    }
    
    # Prepare data.frame for limma (and remove any columns containing 'QC')
    bio_limma_data <- data.frame(assay(input_data)[, !str_detect(colnames(input_data), 'QC')])
    
    # Use shortname as rownames unless rownames provided
    if (is.null(rownames)) {
      rownames(bio_limma_data) <- make.unique(rowData(input_data)$shortname)
    } else {
      rownames(bio_limma_data) <- rownames
    }
    
    # Ensure that NA values are actually NA, and not character type 'NA'
    if (!is.null(metadata_var)) {
      if (!is.null(metadata_condition)) {
        test_df <- data.frame(input_data@metadata$metadata[, metadata_var])[metadata_condition,]
      } else {
        test_df <- data.frame(input_data@metadata$metadata[, metadata_var])
      }
    } else {
      if (!is.null(metadata_condition)) {
        test_df <- data.frame(input_data@metadata$metadata)[metadata_condition,]
      } else {
        test_df <- data.frame(input_data@metadata$metadata)
      }
    }
  }

  # Ensure that character 'NA' values are switched with real 'NA' values
  for (i in 1:ncol(test_df)) {
    if (class(test_df[, i]) == 'character' | class(test_df[, i]) == 'factor') {
      test_df[, i] <- ensure_NA(test_df[, i])
    }
  }
  
  # Attempt to determine which columns are being referenced by the formula
  if (!is.null(model_formula_as_string)) {
    split_formula <- unlist(str_split(model_formula_as_string, pattern = '[+\\d~ ]* '))
    split_formula <- split_formula[!split_formula == ''] # Remove any empty character vector remnants
    
    if (length(split_formula) > 1) {
      test_df <- test_df[, split_formula] # Subset the identified columns from the test data.frame
    } else if (length(split_formula == 1)) {
      test_df <- data.frame(test_df[, split_formula])
      colnames(test_df) <- split_formula
    } else {
      stop('Incorrect dimensions for the data.frame used to generate the model matrix.')
    }
  }
  
  # Reorder factors where requested
  if (!is.null(factor_reorder_list)) {
    for (i in 1:length(factor_reorder_list)) {
      test_df[, names(factor_reorder_list)[i]] <- factor(test_df[, names(factor_reorder_list)[i]],
                                                         levels = factor_reorder_list[[i]])
    }
  }
  
  # Modify continuous variables where requested
  if (!is.null(continuous_modifier_list)) {
    for (i in 1:length(continuous_modifier_list)) {
      test_df[, names(continuous_modifier_list)[i]] <- sapply(test_df[, names(continuous_modifier_list)[i]],
                                                              FUN = continuous_modifier_list[[i]])
    }
  }
  
  # Filter limma data to remove NA values
  bio_limma_data <- bio_limma_data[,!any_NA(test_df)]
  test_df_colnames <- colnames(test_df)
  test_df <- data.frame(test_df[!any_NA(test_df), ])
  colnames(test_df) <- test_df_colnames
  
  # Create design matrix if none provided
  if (is.null(model_matrix)) {
    if (is.null(model_formula_as_string)) {
      if (!is.null(metadata_var)) {
        if (length(metadata_var == 1)) {
          bio_limma_design <- model.matrix(~ 0 + test_df[,1])
          colnames(bio_limma_design) <- levels(test_df[,1])
        } else {
          stop('Please provide a single character value that matches one of the variables in the phyloseq sample data object.')
        }
      } else {
        stop('Please provide either a model matrix or a formula (as a string - using a combination of the metadata variables you selected) for use with limma.
             Alternatively, please provide a single character value to the metadata_var paramter that matches one of the variables in the phyloseq sample data object.')
      }
    } else {
      bio_limma_design <- model.matrix(as.formula(model_formula_as_string), data = test_df)
    }
  } else {
    bio_limma_design <- model_matrix
  }
  
  # Fit the expression matrix to a linear model
  fit <- lmFit(bio_limma_data, bio_limma_design)
  
  # Define a function to handle creation of the contrasts matrix
  make_contrasts_vector <- function(design_matrix) {
    # Define levels and make new lists
    levels <- colnames(design_matrix)
    num_levels <- length(levels)
    ci <- list()
    cj <- list()
    
    # For loops for generate contrast statements
    for (i in 2:num_levels) {
      ci[length(ci) + 1] <- paste0(levels[1], '-', levels[i])
      for (j in (i + 1):num_levels) {
        cj[length(cj) + 1] <- paste0(levels[i], '-', levels[j])
      }
    }
    
    # Unlist elements and remove NA values, then remove the last element (contrasts itself)
    cx <- c(unlist(ci), unlist(cj)); cx <- cx[!str_detect(cx, 'NA')]
    cx <- cx[1:length(cx) - 1]
    
    # Generate contrasts matrix and return it
    cont_matrix <- makeContrasts(contrasts = cx, levels = levels)
    cont_matrix
  }
  
  # Bayes statistics of differential expression
  if(use_contrast_matrix) {
    
    # Create the contrasts matrix if none provided
    if (is.null(contrast_matrix)) {
      cont_matrix <- make_contrasts_vector(bio_limma_design)
    } else {
      cont_matrix <- contrast_matrix
    }
    
    fit2 <- contrasts.fit(fit, contrasts = cont_matrix)
  } else {
    if (is.null(coefficients)) {
      stop('As you have selected to use coefficients instead of a contrast matrix, please assign indices to the 
           coefficients parameter to select coefficients you want to retain for Bayes statistics.')
    } else {
      fit2 <- contrasts.fit(fit, coefficients = coefficients)
    }
  }
  fit2 <- eBayes(fit2, robust = TRUE, trend = FALSE)
  
  # Generate a list of the differentially abundant features
  get_all_topTables <- function(fit2, top = TRUE) {
    # Prepare an empty list object
    all_topTables <- list()
    
    # Run through each of the coefficients and add topTable to the list
    for (coef in 1:dim(fit2$contrasts)[2]) {
      coef_name <- colnames(fit2$contrasts)[coef]
      bio_topTable <- topTable(fit2, number = dim(fit2)[1], adjust.method = adjust_method, coef = coef)
      if (top == TRUE) {
        bio_topTable <- bio_topTable %>%
          filter(adj.P.Val < adj_pval_threshold) %>%
          filter(abs(logFC) >= logFC_threshold)
      }
      all_topTables[[coef_name]] <- bio_topTable
    }
    
    # Return list object
    all_topTables
  }
  
  # Get both significant results and all results lists
  bio_limma_signif <- get_all_topTables(fit2, top = TRUE)
  bio_limma_all <- get_all_topTables(fit2, top = FALSE)
  
  # Generate a venn diagram of the results
  results <- decideTests(fit2, adjust.method = adjust_method, p.value = adj_pval_threshold, lfc = logFC_threshold)
  if (dim(results)[2] <= 4) {
    venn_diagram <- vennDiagram(results)
  } else {
    venn_diagram <- 'No Venn diagram. Greater than 4 comparisons.'
  }
  
  # Add a direction column to the 'bio_limma_all' topTables
  limma_add_categorical_direction <- function(bio_limma_all) {
    for (i in 1:length(bio_limma_all)) {
      bio_limma_all[[i]] <- data.frame(bio_limma_all[[i]]) %>%
        mutate(direction = case_when(
          adj.P.Val < adj_pval_threshold & logFC <= -logFC_threshold ~ 'Decreased',
          adj.P.Val < adj_pval_threshold & logFC >= logFC_threshold ~ 'Increased',
          adj.P.Val >= adj_pval_threshold | abs(logFC) < logFC_threshold ~ 'NS'
        ))
    }
    bio_limma_all
  }
  
  bio_limma_all <- limma_add_categorical_direction(bio_limma_all)
  
  # Assign test_string if provided
  if (!is.null(legend_metadata_string)) {
    test_string <- legend_metadata_string
  } else {
    test_string <- 'Test Variable'
  }
  
  # Create a vector of strings for color legend
  plot_color_names <- function(test_string, bio_limma_all) {
    vector = ''
    if (use_contrast_matrix) {
      for (i in 1:length(bio_limma_all)) {
        contrast_name <- names(bio_limma_all[i])
        contrast_name <- gsub('-', ' vs. ', contrast_name)
        contrast_name <- paste0(test_string, ':\n', contrast_name)
        vector <- c(vector, contrast_name)
      }
    } else {
      for (i in 1:length(bio_limma_all)) {
        contrast_name <- names(bio_limma_all[i])
        contrast_name <- paste0(test_string, ':\n', contrast_name)
        vector <- c(vector, contrast_name)
      }
    }
    vector <- vector[2:length(vector)]
    vector
  }
  
  plot_color_labels <- plot_color_names(test_string, bio_limma_all)
  
  # Define function to plot multiple volcano plots
  limma_volcano_plots <- function(topTable, plot_title, plot_subtitle, plot_xlab, plot_ylab, plot_color_labels) {
    # Create a blank list to hold plots
    plots <- list()
    
    # Make volcano plot for each topTable
    for (i in 1:length(topTable)) {
      plot <- ggplot(topTable[[i]], aes(x = logFC, y = -log10(adj.P.Val))) +
        geom_point(aes(color = direction)) +
        geom_hline(yintercept = -log10(adj_pval_threshold), linetype = 2) +
        geom_vline(xintercept = -logFC_threshold, linetype = 2) +
        geom_vline(xintercept = logFC_threshold, linetype = 2) +
        scale_color_manual(values = c('Decreased' = 'blue', 'Increased' = 'red', 'NS' = 'grey70'), name = plot_color_labels[i]) +
        geom_text_repel(data = topTable[[i]][topTable[[i]]$adj.P.Val < adj_pval_threshold, ],
                        label = rownames(topTable[[i]][topTable[[i]]$adj.P.Val < adj_pval_threshold, ]),
                        size = 2) +
        labs(title = plot_title,
             subtitle = plot_subtitle,
             x = plot_xlab,
             y = plot_ylab)
      
      # Add plot to list
      var_name <- names(topTable[i])
      plots[[var_name]] <- plot
    }
    
    # Return plots
    plots
  }
  
  # Define function to plot bar plots
  limma_bar_plots <- function(topTable, plot_title, plot_subtitle, plot_xlab, plot_color_labels) {
    # Create blank list to hold plots
    plots <- list()
    
    # Make bar plot for each topTable
    for (i in 1:length(topTable)) {
      m <- topTable[[i]] %>%
        mutate(feature = rownames(topTable[[i]])) %>%
        arrange(abs(logFC)) %>%
        slice_tail(n = 40) %>%
        arrange(logFC) %>%
        mutate(suffix = ifelse(str_length(feature) > 50, '...', '')) %>%
        mutate(feature = factor(make.unique(paste0(substr(feature, 1, 50), suffix)), 
                                levels = make.unique(paste0(substr(feature, 1, 50), suffix))),
               direction = ifelse(logFC > 0, 'Increased', 'Decreased')) %>%
        dplyr::select(-suffix)
      
      plot <- ggplot(m, aes(x = logFC, y = feature)) +
        geom_col(aes(fill = direction), alpha = 0.7, width = 0.8) +
        scale_fill_manual(values = c('Decreased' = 'blue', 'Increased' = 'red'), plot_color_labels[i]) +
        labs(title = plot_title,
             subtitle = plot_subtitle,
             x = plot_xlab,
             y = 'Taxa')
      
      # Add plot to list
      var_name <- names(topTable[i])
      plots[[var_name]] <- plot
    }
    
    # Return plots
    return(plots)
  }
  
  # Prepare labelling
  if (is.null(volc_plot_title)) {
    plot_title <- 'Differentially Abundant Features'
  } else {
    plot_title <- volc_plot_title
  }
  
  if (is.null(volc_plot_subtitle)) {
    plot_subtitle <- NULL
  } else {
    plot_subtitle <- volc_plot_subtitle
  }
  
  if (is.null(volc_plot_xlab)) {
    plot_xlab <- 'log2FC'
  } else {
    plot_xlab <- volc_plot_xlab
  }
  
  if (is.null(volc_plot_ylab)) {
    if (adjust_method == 'none') {
      plot_ylab <- 'Significance\n-log10(p-value)'
    } else {
      plot_ylab <- 'Significance\n-log10(adjusted p-value)'
    }
  } else {
    plot_ylab <- volc_plot_ylab
  }
  
  # Define a function to plot significant features individually
  limma_feature_featureplots <- function(topTable) {
    featureplots_list <- list()
    
    # Loop through the items in the test formula
    for (i in names(topTable)) {
      # Retrive the topTable data
      tt <- topTable[[i]]
      
      # Retrieve the input data for just significant features
      bio_sig_init <- bio_limma_data[rownames(bio_limma_data) %in% rownames(tt),] %>%
        rotate_df()
      
      # Check which test formula variable is contained in this part of the topTable, and save name to 'k'
      if (is.null(model_formula_as_string)) {
        k <- metadata_var
      } else {
        k <- split_formula[which(str_detect(i, split_formula))]
      }
      
      # Attach the test variable and rename the test variable column
      if (is.null(model_formula_as_string)) {
        bio_sig <- cbind(bio_sig_init, test_df[, metadata_var])
        colnames(bio_sig) <- c(colnames(bio_sig_init), metadata_var)
      } else {
        bio_sig <- cbind(bio_sig_init, test_df[, k]) 
        colnames(bio_sig) <- c(colnames(bio_sig_init), k)
      }
      
      # Create a blank list to hold these plots
      p_list <- list()
      
      # Determine y-axis label
      box_ylab <- case_when(
        class(input_data) == 'phyloseq' ~ 'Abundance',
        class(input_data) == 'SummarizedExperiment' ~ 'Intensity'
      )
      
      # Loop through the features
      for (fx in rownames(topTable[[i]])) {
        # Create the plots if test variable is of type numeric
        if (class(test_df[,k]) == 'numeric') {
          # Retrieve the adjusted p-value and direction of change
          adj_p_val <- format(round(tt[fx, 'adj.P.Val'], 6), scientific = TRUE)
          log_fc <- round(as.numeric(tt[fx, 'logFC']), 3)
          direction <- ifelse(log_fc > 0, 'red', 'blue')
          
          # Generate subtitle depending on adjust method
          if (adjust_method == 'none') {
            p_subtitle <- paste0('p-value = ', adj_p_val, '; log2FC = ', log_fc)
          } else {
            p_subtitle <- paste0('Adj. p-value = ', adj_p_val, '; log2FC = ', log_fc)
          }
          
          # Generate the plot
          p <- ggplot(bio_sig, aes_string(x = k, y = as.name(fx))) +
            geom_smooth(method = 'lm', formula = y ~ x,se = FALSE, color = direction) +
            stat_poly_eq(method = 'lm', formula = y ~ x, size = 3) +
            geom_point() +
            guides(color = 'none') +
            labs(title = fx,
                 subtitle = p_subtitle,
                 x = str_to_sentence(k),
                 y = box_ylab) +
            theme(text = element_text(size = 8))
          
          p_list[[fx]] <- p
        }
        
        # Create the plots if test variable is of type character or factor
        if (class(test_df[,k]) == 'factor') {
          # Work out the y-positions for manual stats
          p_ydist <- diff(range(bio_sig[, fx])) # Distance between min and max y values
          p_ymax <- max(bio_sig[, fx]) # Max y value
          num_levels <- length(levels(test_df[,k]))
          
          # If using standard Wilcox comparison testing, make the normal plot
          if (redo_boxplot_stats) {
            # Generate comparisons for stat_compare_means
            unique_levels <- unique(levels(test_df[,k]))
            combinations <- rlang::flatten(lapply(seq_along(unique_levels), 
                                           function(x) combn(unique_levels, x, FUN = list)))
            combinations <- combinations[sapply(combinations, length) == 2]
            discard <- str_detect(combinations, 'NA', negate = TRUE)
            combinations <- combinations[discard]
            for (comb in seq_along(combinations)) {
              combinations[[comb]] <- sort(combinations[[comb]])
            }
            
            # Generate the plot
            p <- ggplot(bio_sig, aes_string(x = k, y = as.name(fx))) +
              geom_boxplot(aes_string(fill = k)) +
              geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = p_ydist / 30) +
              stat_compare_means(comparisons = combinations, size = 2.5) +
              scale_fill_jama(alpha = 0.6) +
              guides(fill = 'none') + 
              labs(title = fx,
                   x = str_to_sentence(k),
                   y = box_ylab) +
              coord_cartesian(ylim = c(NA, p_ymax + ((p_ydist / 4) * (num_levels - 1)))) + 
              theme(text = element_text(size = 8))
          } else {
            y_position <- p_ymax + (p_ydist / 20)
            
            # Do stats and feed in values from limma
            stats_manual <- wilcox_test(data = bio_sig, formula = as.formula(paste0('`', fx, '` ~ ', k)), 
                                        ref.group = levels(test_df[,k])[1], 
                                        comparisons = list(c(levels(test_df[,k])[1], gsub(k, '', i)))) %>%
              mutate(limma_padj = format(round(tt[fx, 'adj.P.Val'], 6), scientific = TRUE),
                     y.position = y_position)
            
            # Generate the plot
            p <- ggplot(bio_sig, aes_string(x = k, y = as.name(fx))) +
              geom_boxplot(aes_string(fill = k)) +
              geom_dotplot(binaxis = 'y', stackdir = 'center', binwidth = p_ydist / 30) +
              stat_pvalue_manual(data = stats_manual, label = 'limma_padj', size = 3) +
              scale_fill_jama(alpha = 0.6) +
              guides(fill = 'none') + 
              labs(title = fx,
                   x = str_to_sentence(k),
                   y = box_ylab) +
              coord_cartesian(ylim = c(NA, p_ymax + (p_ydist / 10))) + 
              theme(text = element_text(size = 8))
          }
          p_list[[fx]] <- p
        }
      }
      
      # Add the plots to the feature plots list
      featureplots_list[[i]] <- p_list
    }
    
    # Return the feature plots list
    return(featureplots_list)
  }
  
  # Plot all plot types
  volcano_plots <- limma_volcano_plots(bio_limma_all, plot_title, plot_subtitle, plot_xlab, plot_ylab, plot_color_labels)
  
  bar_plots <- limma_bar_plots(bio_limma_signif, plot_title, plot_subtitle, plot_xlab, plot_color_labels)
  
  feature_plots <- limma_feature_featureplots(bio_limma_signif)
  
  if (!is.null(plot_output_folder)) {
    # Save volcano plots
    for (plot in names(volcano_plots)) {
      if (!is.null(plot_file_prefix)) {
        plot_fp <- here::here(plot_output_folder, paste(plot_file_prefix, plot, 'volcplot.pdf', sep = '_'))
      } else {
        plot_fp <- here::here(plot_output_folder, paste('bio_limma', plot, 'volcplot.pdf', sep = '_'))
      }
      ggsave(plot_fp, volcano_plots[[plot]], width = 20, height = 16, units = 'cm')
    }
    
    # Save bar plots
    for (plot in names(bar_plots)) {
      if (!is.null(plot_file_prefix)) {
        plot_fp <- here::here(plot_output_folder, paste(plot_file_prefix, plot, 'barplot.pdf', sep = '_'))
      } else {
        plot_fp <- here::here(plot_output_folder, paste('bio_limma', plot, 'barplot.pdf', sep = '_'))
      }
      
      if (nrow(bio_limma_signif[[plot]] > 0)) {
        bio_limma_signif_rows <- ifelse(nrow(bio_limma_signif[[plot]]) > 40, 40, nrow(bio_limma_signif[[plot]]))
        barplot_height <- (bio_limma_signif_rows * 0.4) + 2.5
        ggsave(plot_fp, bar_plots[[plot]], width = 20, height = barplot_height, units = 'cm')
      }
    }
    
    # Save feature plots
    for (plot in names(feature_plots)) {
      if (!is.null(plot_file_prefix)) {
        plot_fp <- here::here(plot_output_folder, paste(plot_file_prefix, plot, 'featureplots.pdf', sep = '_'))
      } else {
        plot_fp <- here::here(plot_output_folder, paste('bio_limma', plot, 'featureplots.pdf', sep = '_'))
      }
      
      if (length(feature_plots[[plot]]) > 0) {
        plot_pages <- ggarrange(plotlist = feature_plots[[plot]], nrow = 4, ncol = 3)
        
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
  }
  
  # Get input metadata depending on the input type
  if (class(input_data) == 'phyloseq') {input_metadata <- data.frame(sample_data(input_data))}
  if (class(input_data) == 'SummarizedExperiment') {input_metadata <- data.frame(input_data@metadata$metadata)}
  
  # Make a list of all components above to return from function
  if (use_contrast_matrix) {
    return_list <- list(input_data = bio_limma_data,
                        input_metadata = input_metadata,
                        test_variables = test_df,
                        model_matrix = bio_limma_design,
                        contrast_matrix = cont_matrix,
                        limma_significant = bio_limma_signif,
                        limma_all = bio_limma_all,
                        volcano_plots = volcano_plots,
                        bar_plots = bar_plots,
                        feature_plots = feature_plots,
                        venn_diagram = venn_diagram)
  } else {
    return_list <- list(input_data = bio_limma_data,
                        input_metadata = input_metadata,
                        test_variables = test_df,
                        model_matrix = bio_limma_design,
                        coefficients = coefficients,
                        limma_significant = bio_limma_signif,
                        limma_all = bio_limma_all,
                        volcano_plots = volcano_plots,
                        bar_plots = bar_plots,
                        feature_plots = feature_plots,
                        venn_diagram = venn_diagram)
  }
  
  
  # Return the list
  return(return_list)
}
