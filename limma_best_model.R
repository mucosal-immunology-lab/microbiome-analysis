############################################################################################
# Copyright (c) 2022 - Mucosal Immunology Lab, Monash University, Melbourne, Australia     #
# Author: Matthew Macowan                                                                  #
# This script is provided under the MIT licence (see LICENSE.txt for details)              #
############################################################################################

# Define a function to select the best model formula from a selection of input parameters
limma_best_model <- function(input_data, key_variable = NULL, other_parameters = NULL, selection_metric = 'bic') {
  # Load required packages
  library(pacman)
  pkgs <- c('BiocGenerics', 'base', 'biomformat', 'IRanges', 'S4Vectors', 'dplyr', 'plotly', 
            'ggplot2', 'ggtree', 'ggpubr', 'rstatix', 'Matrix', 'limma', 'phyloseq', 'stats', 
            'stringr', 'utils', 'purrr')
  pacman::p_load(char = pkgs)
  
  # Perform sanity checks
  if (is.null(key_variable)) {stop('Please provide the key variable you are interested in testing downstream.')}
  if (is.null(other_parameters)) {stop('Please provide additional variables for testing for the best model fit.')}
  
  # Define a function to ensure that character 'NA' values are switched with real 'NA' values
  ensure_NA <- function(vector) {
    replace(vector, vector == 'NA', NA)
  }
  
  # Remove samples with NA values for the key_variable
  if (class(input_data) == 'phyloseq') {
    input_data <- prune_samples(!is.na(ensure_NA(sample_data(input_data)[[key_variable]])), input_data)
    
   # Retrieve sample data and OTU table data
    model_matrix_data <- data.frame(sample_data(input_data))
    count_data <- data.frame(otu_table(input_data))
  }
  if (class(input_data) == 'SummarizedExperiment') {
    input_data <- input_data[, !is.na(ensure_NA(input_data@metadata$metadata[[key_variable]]))]
    input_metadata <- input_data@metadata$metadata[!is.na(ensure_NA(input_data@metadata$metadata[[key_variable]])),]
    
    if (class(input_metadata[[key_variable]]) == 'factor') {
      input_metadata[[key_variable]] <- droplevels(input_metadata[[key_variable]])
    }
    
    # Retrieve sample data and counts data
    model_matrix_data <- data.frame(input_metadata)
    count_data <- data.frame(assay(input_data))
  }
  
  # Prepare the combinations of test parameters provided
  combinations <- flatten(lapply(seq_along(c(key_variable, other_parameters)), function(x) combn(c(key_variable, other_parameters), x, FUN = list)))
  
  for (i in seq_along(combinations)) {
    combinations[[i]] <- paste(combinations[[i]], collapse = ' + ')
  }
  
  combinations <- unlist(combinations)
  
  # Prepare model matrixes for each combination
  model_list <- list()
  for (i in seq_along(combinations)) {
    formula <- as.formula(paste0('~', combinations[i]))
    model <- model.matrix(formula, data = model_matrix_data)
    model_list[[i]] <- model
  }
  
  # Check for the best model
  model_scores <- data.frame(table(limma::selectModel(count_data, model_list, criterion = selection_metric)$pref),
                             model = combinations) %>%
    mutate(best = ifelse(Freq == max(Freq), 'Best', 'Other')) %>%
    arrange(desc(Freq)) %>%
    mutate(formula_length = unlist(lapply((str_split(as.character(model), ' \\+ ')), length))) %>%
    mutate(model = factor(model, levels = model),
           freq_best_perc = Freq / sum(Freq) * 100)
  colnames(model_scores) <- c('model_number', 'freq_as_best', 'Model', 'Best', 'formula_length', 'freq_best_perc')
  
  model_scores2 <- model_scores %>%
    arrange(freq_best_perc) %>%
    mutate(Model = factor(Model, levels = Model))
  
  if (nrow(model_scores2) > 10) {
    model_scores2 <- slice_tail(model_scores2, n = 10)
  }
  
  plot_subtitle <- case_when(
    selection_metric == 'aic' ~ 'Akaike\'s Information Criterion (AIC)',
    selection_metric == 'bic' ~ 'Bayesian Information Criterion (BIC)'
  )
  
  # Plot the AIC scores
  scores_plot <- ggplot(model_scores2, aes(x = freq_best_perc, y = Model)) +
    geom_col(aes(fill = Best)) +
    scale_fill_manual(values = c('Best' = '#1261A0', 'Other' = 'grey60')) +
    guides(fill = 'none') +
    labs(title = 'Best Fitting Model',
         subtitle = plot_subtitle,
         x = 'Model is Best for this\nPercentage of Features')
  print(scores_plot)
  
  # Retrieve the best model
  best_model <- paste0('~ ', model_scores$Model[1])
  if (!str_detect(best_model, key_variable)) {
    best_model <- paste0(best_model, ' + ', key_variable)
  }
  best_model_matrix <- model.matrix(as.formula(best_model), data = model_matrix_data)
  
  # Return data or graph
  return(list(formula_string = best_model,
              key_var_coef = which(str_detect(colnames(best_model_matrix), key_variable)),
              scores_plot = scores_plot,
              model_scores = model_scores,
              selection_metric = str_to_upper(selection_metric)))
}