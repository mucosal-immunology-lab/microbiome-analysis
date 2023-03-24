############################################################################################
# Copyright (c) 2022 - Respiratory Immunology Lab, Monash University, Melbourne, Australia #
# Author: Matthew Macowan                                                                  #
# This script is provided under the MIT licence (see LICENSE.txt for details)              #
############################################################################################

# Define a function to process the raw output from kraken2 and prepare it for use in a phyloseq object
kraken2_preprocess <- function(filepath = NULL, otu_id_colname = '#OTU ID', lineage_colname = 'Consensus Lineage') {
  # Load required packages
  pkgs <- c('data.table', 'tidyverse')
  for (pkg in pkgs) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
  
  # Import the raw bracken abundance data using the provided file path
  kraken2_input <- read_csv(file = filepath, skip = 1)
  
  ### PROCESS DATA FOR THE PHYLOSEQ TAXONOMY TABLE ###
  # Extract the consensus lineage information from the raw data
  consensus_lineage_kraken2 <- data.frame('otu_id' = kraken2_input[[otu_id_colname]],
                                          'lineage' = kraken2_input[[lineage_colname]])
  
  # Extract lineage information into column
  consensus_lineage_kraken2 <- consensus_lineage_kraken2 %>%
    mutate(kingdom = gsub('k__([a-zA-Z ]*).*', '\\1', consensus_lineage_kraken2$lineage),
           phylum = gsub('.*p__([a-zA-Z ]*).*', '\\1', consensus_lineage_kraken2$lineage),
           class = gsub('.*c__([a-zA-Z ]*).*', '\\1', consensus_lineage_kraken2$lineage),
           order = gsub('.*o__([a-zA-Z ]*).*', '\\1', consensus_lineage_kraken2$lineage),
           family = gsub('.*f__([a-zA-Z ]*).*', '\\1', consensus_lineage_kraken2$lineage),
           genus = gsub('.*g__([a-zA-Z ]*).*', '\\1', consensus_lineage_kraken2$lineage),
           species = gsub('.*s__([a-zA-Z ]*).*', '\\1', consensus_lineage_kraken2$lineage)) %>%
    dplyr::select(-lineage)
  consensus_lineage_kraken2$species <- replace(consensus_lineage_kraken2$species, consensus_lineage_kraken2$species == 'sp', NA)
  consensus_lineage_kraken2$species <- replace(consensus_lineage_kraken2$species, consensus_lineage_kraken2$species == '', NA)
  
  # Change unassigned taxonomic classifications to 'Unknown'
  consensus_lineage_kraken2$phylum <- replace(consensus_lineage_kraken2$phylum, consensus_lineage_kraken2$phylum == '', 'Unknown')
  consensus_lineage_kraken2$class <- replace(consensus_lineage_kraken2$class, consensus_lineage_kraken2$class == '', 'Unknown')
  consensus_lineage_kraken2$order <- replace(consensus_lineage_kraken2$order, consensus_lineage_kraken2$order == '', 'Unknown')
  consensus_lineage_kraken2$family <- replace(consensus_lineage_kraken2$family, consensus_lineage_kraken2$family == '', 'Unknown')
  consensus_lineage_kraken2$genus <- replace(consensus_lineage_kraken2$genus, consensus_lineage_kraken2$genus == '', 'Unknown')
  
  # Give each taxa a unique ID based on OTU ID
  consensus_lineage_kraken2$taxa <- paste(consensus_lineage_kraken2$genus, consensus_lineage_kraken2$species, consensus_lineage_kraken2$otu_id, sep = '_')
  
  # Move OTU IDs to rownames
  consensus_lineage_kraken2 <- consensus_lineage_kraken2 %>%
    column_to_rownames(var = 'otu_id')
  
  ### PROCESS DATA FOR THE PHYLOSEQ OTU TABLE ###
  # Use OTU ID temporarily for taxon rownames, and extract only read counts
  kraken2_raw <- kraken2_input %>%
    column_to_rownames(var = otu_id_colname)
  kraken2_raw[[lineage_colname]] <- NULL
  
  ### COMBINE PREPARED ELEMENTS INTO A LIST FOR RETURN FROM THE FUNCTION ###
  return_list <- list(kraken2_tax_table = consensus_lineage_kraken2,
                      kraken2_otu_table = kraken2_raw)
  
  return(return_list)
}
