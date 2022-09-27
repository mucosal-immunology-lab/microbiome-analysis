# Function to prepare and save OTU tables
otu_to_csv <- function(phyloseq_object, output_filename) {
  otu_table <- data.frame(otu_table(phyloseq_object))
  tax_table <- data.frame(tax_table(phyloseq_object))
  colnames(otu_table) <- rownames(sample_data(phyloseq_object))
  otu_table <- rownames_to_column(otu_table, var = 'sequence')
  otu_table <- cbind(tax_table, otu_table)
  write_csv(otu_table, output_filename)
}