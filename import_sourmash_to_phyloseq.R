library(tidyverse)
library(phyloseq)

#load metadata and reformat
metadata <- read_excel('davis_water_metadata.xlsx') %>% 
  select(c(renamed,treatment,temperature)) %>%
  rename(sample = renamed)%>%
  column_to_rownames('sample')
# read in the sourmash taxonomy results from all samples into a single dataframe
sourmash_tax_results <- Sys.glob("sourmash_gather_with_lineages/*.with-lineages.csv") %>%
  map_dfr(read_csv, col_types = "ddddddddcccddddcccdc") %>%
  mutate(name = gsub(" .*", "", name))

# We need two tables: a tax table and an "otu" table. 
# The tax table will hold the taxonomic lineages of each of our gather matches.
# To make this, we'll make a table with two columns: the genome match and the lineage of the genome.
# The "otu" table will have the counts of each genome in each sample.
# We'll call this our gather_table.

tax_table <- sourmash_tax_results %>%
  select(name, lineage) %>%
  distinct() %>%
  separate(lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  column_to_rownames("name")


gather_table <- sourmash_tax_results %>% 
  mutate(n_unique_kmers = (unique_intersect_bp / scaled) * average_abund) %>% # calculate the number of uniquely matched k-mers and multiply it by average k-mer abundance
  select(query_name, name, n_unique_kmers) %>% # select only the columns that have information we need
  pivot_wider(id_cols = name, names_from = query_name, values_from = n_unique_kmers) %>% # transform to wide format
  replace(is.na(.), 0) %>% # replace all NAs with 0 
  column_to_rownames("name") # move the metagenome sample name to a rowname

# create phyloseq object
physeq_all_reads <- phyloseq(otu_table(gather_table, taxa_are_rows = T),
                             tax_table(as.matrix(tax_table)),
                             sample_data(metadata))
saveRDS(physeq_all_reads,'physeq_all_reads.rds')

#------------------------------------------------------------------------------------------------------
#now we will perform the same thing for arg containing contigs

sourmash_arg_contigs_tax_results <- Sys.glob("arg_contig_sourmash_gather_with_lineage/*.with-lineages.csv") %>%
  map_dfr(read_csv, col_types = "ddddddddcccddddcccdc") %>%
  mutate(name = gsub(" .*", "", name))
arg_tax_table <- sourmash_arg_contigs_tax_results %>%
  select(name, lineage) %>%
  distinct() %>%
  separate(lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
  column_to_rownames("name")
arg_gather_table <- sourmash_arg_contigs_tax_results %>% 
  mutate(n_unique_kmers = (unique_intersect_bp / scaled) * average_abund) %>% 
  select(query_name, name, n_unique_kmers) %>% 
  pivot_wider(id_cols = name, names_from = query_name, values_from = n_unique_kmers) %>% 
  replace(is.na(.), 0) %>% 
  column_to_rownames("name")
physeq_arg_contigs <- phyloseq(otu_table(arg_gather_table, taxa_are_rows = T),
                               tax_table(as.matrix(arg_tax_table)),
                               sample_data(metadata))
saveRDS(physeq_arg_contigs,'physeq_arg_contigs.rds')

