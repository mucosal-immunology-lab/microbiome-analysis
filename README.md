# Microbiome Analysis

Here we will provide a selection of analytical tools and utilities for the processing of microbiome data derived from either 16S rRNA amplicon sequencing or shotgun metagenomics.

## Table of Contents

  - [Processing of raw data](#processing-of-raw-data)
  - [Downstream preprocessing - 16S rRNA amplicon sequencing](#downstream-preprocessing---16s-rrna-amplicon-sequencing)
    - [Load bacterial 16S data](#load-bacterial-16s-data)
    - [Prepare metadata](#prepare-metadata)
    - [Import data into a phyloseq object](#import-data-into-a-phyloseq-object)
  - [Downstream preprocessing - Shotgun metagenomics](#downstream-preprocessing---shotgun-metagenomics)
    - [Load Kraken2/Bracken output](#load-kraken2bracken-output)
    - [Prepare metadata](#prepare-metadata-1)
    - [Import data into a phyloseq object](#import-data-into-a-phyloseq-object-1)
  - [Filtering and normalisation](#filtering-and-normalisation)
    - [Diversity](#diversity)
    - [Normalisation](#normalisation)
    - [Write OTU tables to `.csv` files](#write-otu-tables-to-csv-files)
    - [Data agglomeration](#data-agglomeration)
  - [Alpha diversity](#alpha-diversity)
    - [Example: longitudinal experiment with two groups](#example-longitudinal-experiment-with-two-groups)
  - [Beta diversity](#beta-diversity)

## Processing of raw data

Please refer to the relevant pipeline for processing of raw sequencing reads:

* 16S rRNA amplicon sequencing: [DADA2 pipeline](https://github.com/respiratory-immunology-lab/microbiome-dada2)
* Shotgun metagenomics: [Sunbeam pipeline](https://github.com/respiratory-immunology-lab/microbiome-shotgun)

## Downstream preprocessing - 16S rRNA amplicon sequencing

After completing the steps in the DADA2 pipeline, you should have three input files for generation of a `phyloseq` object (a container object to hold your taxa counts, information, and sample data):

* `seqtab_nochim.rds`: the counts data for each ASV after the removal of chimeras.
* `taxonomy_species.rds`: the taxonomic lineage information for each ASV, as assigned using the SILVA database.
* `tree.rds`: a *de novo* phylogenetic tree generated in the DADA2 pipeline that encompasses all of your ASVs.

### Load bacterial 16S data

Firstly, load in your input files. You may wish to change the sample names (i.e. column names) of your `seqtab_nochim.rds` file, but be sure that the sample names for your metadata at the next step match (or are changed to match).

The `highest_ClassID()` function is provided as a script [here](./highest_ClassID.R)

```r
# Load bacterial count table seqtab_nochim.rds and corresponding taxonomic classification taxonomy_species.rds
bact_seqtab <- t(data.frame(readRDS(here::here('data', 'dada2', 'seqtab_nochim.rds'))))

bact_taxonomy <- data.frame(readRDS(here::here('data', 'dada2', 'taxonomy_species.rds')))
bact_tree <- readRDS(here::here('data', 'dada2', 'tree.rds'))

# Add a taxonomy column 'ID' to the taxonomy for visualisation purposes (highest level of classification and a unique ASV number)
bact_taxonomy$ID <- highest_ClassID(bact_taxonomy)
```

### Prepare metadata

Your metadata should be prepared with columns in the same order as your OTU table, with identical column names.

To ensure everything is in the correct order, we can also match the sequencing data with the metadata and subset the dataset to keep only `bact_seqtab` samples with corresponding metadata.

```r
# Keep only metadata rows that match samples in the bact_seqtab
metadata <- metadata[which(rownames(metadata) %in% colnames(bact_seqtab)),]

# Now, keep only bact_seqtab samples with corresponding metadata
bact_seqtab <- bact_seqtab[, which(colnames(bact_seqtab) %in% rownames(metadata))]

# Take a look at the dimensions of the remaining seqtab
dim(bact_seqtab)
```

### Import data into a phyloseq object

`phyloseq` is an R package to import, store, analyse, and graphically display complex phylogenetic sequencing data that has already been processed.

An example for its creation is given below.

```r
# Get data ready for phyloseq creation
otu_mat <- as.matrix(bact_seqtab[, sort(colnames(bact_seqtab))])
tax_mat <- as.matrix(bact_taxonomy)
samples_df <- metadata[colnames(bact_seqtab),]

# Use data to create the individual phyloseq object elements
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(samples_df)
sample_names(samples) <- colnames(otu_mat)

# Import data into a phyloseq object
bact_data_raw <- phyloseq(OTU, TAX, samples)

# Remove temp files
rm(list = c('otu_mat', 'tax_mat', 'samples_df', 'OTU', 'TAX'))
```

## Downstream preprocessing - Shotgun metagenomics

Given that much of the quality control and filtering of our shotgun metagenomic data has already taken place using the Sunbeam pipeline, and we also have taxonomic assignments using the combined power of the Kraken2 and Bracken tools, we do not require as much local pre-processing of data in R before downstream analysis.

### Load Kraken2/Bracken output

The first step is to load in your `all_samples_kraken2.csv` or `all_samples_bracken.csv` file, and extract the OTU ID and consensus lineage information. This will help set you up for preparation of a `phyloseq` container object to hold your data.

We provide a custom function here ([`kraken2_preprocess.R`](./kraken2_preprocess.R) for producing a taxonomy table and OTU table with unique identifiers from your output file. You may wish to correct the column names after this function however.

The `kraken2_preprocess()` function will return a list with two elements: firstly an object called `kraken2_tax_table` which holds the taxonomy table, and secondly one called `kraken2_otu_table` which holds the OTU read count information (these can be individually extracted using the typical `$` notation).

```r
# Read in and pre-process raw kraken2 data
bact_kraken2_input <- kraken2_preprocess(filepath = here::here('input', 'shotgun_data', 'all_samples_kraken2.csv'))

# Extract individual elements
kraken2_otu_table <- bact_kraken2_input$kraken2_otu_table
kraken2_tax_table <- bact_kraken2_input$kraken2_tax_table
```

### Prepare metadata

Your metadata should be prepared with columns in the same order as your OTU table, with identical column names.

### Import data into a phyloseq object

`phyloseq` is an R package to import, store, analyse, and graphically display complex phylogenetic sequencing data that has already been processed.

An example for its creation is given below.

```r
# Get data ready for phyloseq creation
otu_mat <- as.matrix(kraken2_otu_table)
tax_mat <- as.matrix(kraken2_tax_table)
samples_df <- sample_metadata

# Use data to create the individual phyloseq object elements
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(samples_df)
sample_names(samples) <- colnames(otu_mat)

# Import data into a phyloseq object
bact_data_raw <- phyloseq(OTU, TAX, samples)

# Retain only bacterial taxa (you could also choose another kingdom here)
bact_data_raw <- subset_taxa(bact_data_raw, kingdom == 'Bacteria')

# Remove temp files
rm(list = c('otu_mat', 'tax_mat', 'samples_df', 'OTU', 'TAX'))
```

## Filtering and normalisation

At this stage, we should now have a single, unified `phyloseq` object that contains all of our count data, taxonomic information, and sample metadata:

* `bact_data_raw`: the raw `phyloseq` object before filtering and normalisation (can be derived from either 16S rRNA amplicon or shotgun metagenomic sequencing data).

We can now perform some filtration steps to remove samples with low read counts, and also taxa with very few reads or those that are only present in a small number of samples. 

```r
# Remove samples with less than [minreadsThreshold] reads
minreadsThreshold <- 10000 # originally set to 10000
bact_data_samples <- prune_samples(colSums(otu_table(bact_data_raw)) > minreadsThreshold, bact_data_raw)

# Filter out taxa based on a minimum number of reads [detectionThreshold] and prevalence [prevalenceThreshold]
detectionThreshold <- 0 # originally set to 0 (i.e. 0 reads minimum)
prevalenceThreshold <- 0.1 # originally set to 0.1 (i.e. found in 10% of samples)
bact_data_filtered <- core(bact_data_samples,
                           detection = detectionThreshold,
                           prevalence = prevalenceThreshold,
                           include.lowest = FALSE)

# Remove samples with less than [minreadsThreshold] reads, after the prevalence filtering.
minreadsThreshold <- 5000 # originally set to 5000
bact_data_filtered <- prune_samples(colSums(otu_table(bact_data_filtered)) > minreadsThreshold, bact_data_filtered)

# Generate a summary of filtration
data_summary <- data.frame('Summary' = c('OTUs/genera', 'Samples'),
                           'Bacteria_before_filtering' = dim(otu_table(bact_data_raw)),
                           'Bacteria_after_filtering' = dim(otu_table(bact_data_filtered)))

# Print the table as a kable object
kable(data_summary) %>%
  kable_styling(bootstrap_options = 'striped', full_width = FALSE)

# Save the filtered phyloseq to an .rds file
saveRDS(bact_data_filtered, here::here('output', 'bact_data_filtered.rds'))
```

### Diversity

We can add Shannon diversity index information at this time, so it will be incorporated into any normalised datasets.

```{r}
# Estimate Shannon index
sample_data(bact_data_filtered)$Diversity <- estimate_richness(bact_data_filtered, split = TRUE, measures = c('Shannon'))$Shannon
```

### Normalisation

Prior to conducting multivariate analysis, it is important to consider the structure of the data. First, high-throughput sequencing data has high variability in library size and therefore differences between samples do not demonstrate true biological diversity. Second, shotgun counts are sparse, with most components containing a zero count. To accommodate for the high variability and sparseness, we require normalisation techniques to improve downstream statistical analysis.

We will use 3 different normalisation methods:

**Centered log ration transformation (CLR)**: 
This normalisation is robust to compositionality. This is the preferred transformation when using some multivariate approaches such as sPLS-DA.

**Cumulative Sum Scaling transformation (CSS)**: 
This normalisation is robust to compositionality and has been specifically developed for microbiome data.

**Log Cumulative Sum Scaling transformation (logCSS)**: 
The log + 1 transformation of CSS-normalised data.

```r
# Perform centred log ratio transformation
bact_data_CLR <- microbiome::transform(bact_data_filtered, transform = 'clr')
saveRDS(bact_data_CLR, here('output', 'bact_data_CLR.rds'))

# Perform cumulative sum scaling transformation
bact_data_CSS <- bact_data_filtered
otu_table(bact_data_CSS) <- 
  otu_table(MRcounts(cumNorm(phyloseq_to_metagenomeSeq(bact_data_filtered), p = 0.5)), taxa_are_rows = TRUE)
saveRDS(bact_data_CSS, here('output', 'bact_data_CSS.rds'))

# Perform cumulative sum scaling tranformation followed by log transformation
bact_data_logCSS <- microbiome::transform(bact_data_CSS, transform = 'log')
saveRDS(bact_data_logCSS, here('output', 'bact_data_logCSS.rds'))
```

### Write OTU tables to `.csv` files

At this point, we can also write each of the OTU tables to `.csv` files using the `otu_to_csv()` function provided [here](./otu_to_csv.R).

```r
# Write OTU tables to .csv file
otu_to_csv(bact_data_filtered, here::here('output', 'otu_table_bact_data_filtered.csv'))
otu_to_csv(bact_data_CLR, here::here('output', 'otu_table_bact_data_CLR.csv'))
otu_to_csv(bact_data_CSS, here::here('output', 'otu_table_bact_data_CSS.csv'))
otu_to_csv(bact_data_logCSS, here::here('output', 'otu_table_bact_data_logCSS.csv'))
```

### Data agglomeration

We may also want to generate agglomerated forms of the `phyloseq` at different taxonomic levels. This will allow us to interrogate differences at the genus or family levels, for example.

*Double check whether your taxonomy table uses the capitalised form of the taxonomic level, i.e. 'genus' or 'Genus'.*

```r
# Prepare genus-agglomerated phyloseq objects
bact_genus_filtered <- tax_glom(bact_data_filtered, taxrank = 'Genus')
bact_genus_logCSS <- tax_glom(bact_data_logCSS, taxrank = 'Genus')

# Save objects to RDS files
saveRDS(bact_genus_filtered, here::here('output', 'bact_genus_filtered.rds'))
saveRDS(bact_genus_logCSS, here::here('output', 'bact_genus_logCSS.rds'))

# Write OTU tables to .csv file
otu_to_csv(bact_genus_filtered, here::here('output', 'otu_tables', 'otu_table_bact_genus_filtered.csv'))
otu_to_csv(bact_genus_logCSS, here::here('output', 'otu_tables', 'otu_table_bact_genus_logCSS.csv'))
```

## Alpha diversity

Alpha diversity represents **within-sample** diversity. We can now make use of the Shannon diversity metric information we added to our phyloseq objects earlier to investigate whether there are changes based on some sample metadata of our choosing. 

The approach for this will be different depending on whether your metadata of interest is categorical or numerical, and how you want to split up your data (if at all &ndash; for example, you may have a numeric explanatory variable, but also a categorical grouping factor).

### Example: longitudinal experiment with two groups

Perhaps we have a longitudinal study where faecal samples have been taken from the same animals at five distinct time points, and the animals belong to one of two treatment groups, e.g. a Sham group and a drug treatment group.

Firstly we may want to determine if there are changes from the baseline alpha diversity (timepoint `0`). One way we could investigate that is as follows:

```r
# Transform the sample data into a data.frame
bact_data_samples <- data.frame(sample_data(bact_data_filtered))

# Plot the diversity over time by group
(bact_diversity_plot <- ggplot(bact_data_samples, aes(x = days_postTx, y = Diversity)) +
    geom_path(aes(group = rat_id), col = 'grey70') + # draw a path between sample from the same animal
    geom_violin(aes(y = Diversity, fill = group), scale = 'area') +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, binwidth = 1/30 * diff(range(bact_data_samples$Diversity))) +
    geom_boxplot(aes(y = Diversity, fill = group), width = 0.2) +
    scale_fill_jama(alpha = 0.5, name = 'Treatment Group') +
    guides(fill = 'none') +
    labs(title = 'Bacterial Diversity (Shannon Index)',
         x = 'Time Since Treatment (days)') +
    facet_grid(cols = vars(group)) +
    stat_compare_means(comparisons = list(c('0', '7'),
                                          c('0', '14'),
                                          c('0', '21'),
                                          c('0', '28')))
)
ggsave(here::here('figures', 'diversity', 'bact_diversity_by_group.pdf'), bact_diversity_plot,
       width = 16, height = 12, units = 'cm')
```

That would produce this plot:

<p align="center">
    <img src="./assets/bact_diversity_by_group.jpg" width=75%>
</p>

Alternatively, we may decide we want to look at whether the treatment affects alpha diversity between groups over time. For example, there may be some temporal element affecting alpha diversity that is independent of the treatment. This is easy to alter by changing the faceting variable, and removing the `geom_path()` element.

```r
# Transform the sample data into a data.frame
bact_data_samples <- data.frame(sample_data(bact_data_filtered))

# Prepare facet labels (more descriptive than the numbers alone)
# We will use the labeller within facet_grid to use these labels
facet_col_labels <- c('0' = '0 days post-Tx',
                      '7' = '7 days post-Tx',
                      '14' = '14 days post-Tx',
                      '21' = '21 days post-Tx',
                      '28' = '28 days post-Tx')

# Plot the diversity over time between groups
(bact_diversity_plot2 <- ggplot(bact_data_samples, aes(x = group, y = Diversity)) +
    geom_violin(aes(y = Diversity, fill = group), scale = 'area') +
    geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, 
                 binwidth = 1/30 * diff(range(bact_data_samples$Diversity))) +
    geom_boxplot(aes(y = Diversity, fill = group), width = 0.2) +
    scale_fill_jama(alpha = 0.5, name = 'Treatment Group') +
    guides(fill = 'none') +
    labs(title = 'Bacterial Diversity (Shannon Index)',
         x = 'Treatment Group') +
    facet_grid(cols = vars(days_postTx), scales = 'free', labeller = as_labeller(facet_col_labels)) +
    stat_compare_means(comparisons = list(c('Sham', 'DrugXYZ')))
)
ggsave(here::here('figures', 'diversity', 'bact_diversity_by_time.pdf'), bact_diversity_plot2,
       width = 16, height = 12, units = 'cm')
```

That would produce this plot:

<p align="center">
    <img src="./assets/bact_diversity_by_time.jpg" width=75%>
</p>

## Beta diversity

Beta diversity represents **between-sample** diversity. Continuing on from our example above, because the taxonomic composition may vary greatly with time, if we were to plot ordination plots for all samples together, we may lose clarity and minimise resolution/separation of data points at each of the time points individually.

In this example, the data came from 16S rRNA amplicon sequencing, so we have access to the *de novo* phylogenetic tree, and can use the UniFrac distance metric. UniFrac distances make use of the tree, by assigning short distances to a pair of ASVs that are highly similar, and larger distances to those that are more dissimilar.
If however we were using a `phyloseq` derived from shotgun metagenomic sequencing data, we could use another distance metric, e.g. Bray-Curtis (`'bray'`).

To save time and code, let's loop through the ordinations at each timepoint by creating temporary `phyloseq` subsets, add their plotted ordinations to a plot list, and finally combine the plots into a single figure with `ggarrange()`.

*If you just have a single categorical variable however, you can make use of the plot code within the loop, i.e. `p <- plot_ordination()` etc.*

```r
# Choose parameters for downstream analyses
met <- 'PCoA'
dist <- 'unifrac'
bact <- bact_data_logCSS

# Bacterial ordination (separately for each temporal group)
ord_list <- list() # blank list to hold ordination plots
for (day in levels(bact@sam_data$days_postTx)) { # loop through the unique values of 'days_postTx' in the phyloseq sample_data()
  bact_tmp <- prune_samples(bact@sam_data$days_postTx == day, bact) # create temporary subsets of the phyloseq object
  
  p <- plot_ordination(bact_tmp, ordinate(bact_tmp, met, dist, weighted = TRUE), # ordinate via the phyloseq plot_ordination() function
                       title = paste0(day, ' days post-treatment')) +
    stat_ellipse(aes(fill = group), geom = 'polygon', type = 't', level = 0.95, alpha = 0.2) + # add group ellipses
    scale_shape_identity() +
    geom_point(aes(fill = group), shape = 21, size = 3) + # add the individual data points
    scale_fill_jama(name = 'Treatment Group') # scale the fill colour using ggsci scale_fill_jama() function
  
  ord_list[[day]] <- p # add the plot to the plot list
}

# Arrange ordination plots using ggpubr ggarrange() function
ord_plots <- ggarrange(plotlist = ord_list, nrow = 2, ncol = 3, common.legend = TRUE)
(ord_plots <- annotate_figure(ord_plots,
                              top = text_grob(label = paste0('Bacteria ', met, ' ', dist, ' ordination'))))
ggsave(here::here('figures', 'ordination', 'bact_PCoA_ordination_by_time.pdf'), ord_plots, # Save combined figure
       width = 20, height = 16, units = 'cm')
```

The resulting output file looks like this:

<div align="center">
  <img src="./assets/bact_PCoA_ordination_by_time.jpeg" width = 75%>
</div>

## Taxonomic composition

To investigate sample taxonomic composition, we can either look at the average composition of groups (using some sample metadata variable) or at the composition of each individual sample.

In either case, we will first need to filter the taxa using the `microbiome::core()` function to increase the clarity of the resulting barplots. We will place that within the `phyloseq::tax_glom()` function &ndash; to agglomerate to the genus level for example &ndash; and then transform the data to compositional form, and use the `phyloseq::psmelt()` function to pivot the data longer for plotting with `ggplot2`.

### Data preparation

```r
# Filter taxa to increase clarity of barplots
bact_data_tax <- psmelt(microbiome::transform(tax_glom(core(bact_data_filtered, detection = 10, prevalence = 0.2), 
                                                       'Genus', NArm = T), transform = 'compositional', target = 'OTU'))
```

Next, we will prepare a colour vector for the barplots using the `RColorBrewer` package.

```r
# Prepare colours for plotting
library(RColorBrewer)
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
```

### Taxonomic barplot by group

We will first plot the "average" barplot using a grouping variable. This could be a treatment group, time frame, etc. In this example, we will continue on with the drug treatment example from above, using the treatment as our grouping variable and also facet by time frame (by passing the `facet_col_labels` vector we created above to the `labeller` argument of `facet_grid()`).

```r
# Plot taxonomic composition by group over time
(bact_taxonomy_gen <- ggplot(bact_data_tax, aes(x = group, y = Abundance, fill = Genus, NArm = FALSE)) +
    geom_bar(stat = 'identity', position = 'fill', size = 0.01) +
    scale_fill_manual(values = col_vector) + # use the colour vector we created
    labs(title = 'Genera Distribution vs Treatment Group',
         x = 'Treatment Group',
         y = 'Abundance') +
    theme(legend.key.size = unit(3.5, 'mm'), # change the size of each legend box
          legend.position = 'bottom',
          plot.margin = margin(0.25, 3, 0.25, 1.5, unit = 'cm'), # increase margins around the plot
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) + # rotate x-axis labels
    facet_grid(cols = vars(days_postTx), labeller = as_labeller(facet_col_labels))
)
ggsave(here::here('figures', 'composition', 'bact_taxonomy_genera_by_time.pdf'), bact_taxonomy_gen,
       width = 20, height = 15, units = 'cm')
```

The resulting output file looks like this:

<div align="center">
  <img src="./assets/bact_taxonomy_genera_by_time.jpg" width = 75%>
</div>

### Taxonomic composition by individual

Seeing as we also know the IDs of the animals, we can also plot the taxonomic composition for each individual over time, facetting by both group and time frame.

*We also need to add the treatment names to the end of the `facet_col_labels` vector here, because we are now facetting by two factors. The `labeller` just takes a single named vector for all labels.*

```r
# Barplot of abundance for each individual by group over time
(bact_taxonomy_gen_rat <- ggplot(bact_data_tax, aes(x = rat_id, y = Abundance, fill = Genus, NArm = FALSE)) +
    geom_bar(stat = 'identity', position = 'fill', size = 0.01) +
    scale_fill_manual(values = col_vector) +
    labs(title = 'Genera Distribution vs Treatment Group',
         x = 'Rat',
         y = 'Abundance') +
    theme(legend.key.size = unit(3.5, 'mm'),
          legend.position = 'bottom',
          plot.margin = margin(0.25, 3, 0.25, 2.25, unit = 'cm'),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    facet_grid(row = vars(days_postTx), cols = vars(group), scales = 'free', labeller = as_labeller(c(facet_col_labels,
                                                                                                      'Sham' = 'Sham',
                                                                                                      'DrugXYZ' = 'DrugXYZ')))
)
ggsave(here('figures', 'composition', 'bact_taxonomy_gen_rat_by_time.pdf'), bact_taxonomy_gen_rat,
       width = 20, height = 20, units = 'cm')
```

The resulting output file looks like this:

<div align="center">
  <img src="./assets/bact_taxonomy_gen_rat_by_time.jpg" width = 75%>
</div>