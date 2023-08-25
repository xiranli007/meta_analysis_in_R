library(tidyverse)
library(readxl)
library(openxlsx)
library(stringr)

# load AMR_analytic_matrix genertaed from resistom analyzer
amr_hierarchy <- c("header", "type", "class", "mechnism",'group', NA)
amr_matrix <- read_excel('AMR_analytic_matrix.xlsx', sheet = "sorted") %>% 
  pivot_longer(-sample, names_to='ARG', values_to = 'hits') %>%
  separate(ARG, amr_hierarchy,"\\|")
# load metadata
metadata <- read_excel('davis_water_metadata.xlsx')
# load metaxt summary flies
file.names <- dir('metaxa_summaries/', pattern = ".summary.txt", recursive = TRUE, full.names = TRUE)
# a function to extract number of 16S rRNA from metaxa summary files
extract_info <- function (file){
  name <-basename(file) %>% gsub("_metaxa_out.summary.txt", "",., fixed = TRUE)
  num_of_bacterial_rRNA <- read_lines(file, skip=25, n_max=1) %>% parse_number()
  return (c(name, num_of_bacterial_rRNA))
}
#read the metaxa summary files for each sample and extract the information
info <- lapply(file.names,extract_info)
# create metaxa df and organize it
metaxa_df <- as.data.frame(do.call(rbind,info))
colnames(metaxa_df) <- c ('sample','num_of_16srRNA')
metaxa_df$num_of_16srRNA <- as.numeric (metaxa_df$num_of_16srRNA)
# load gene length frome megares library
gene_length_df <- read_excel ('megares_database_v3.00.gene_length.xlsx', col_names = FALSE)
colnames(gene_length_df ) <- c('header',"gene_length")
# combine metaxa, amr_analtic_matrix and metadata,gene_length_df to single big df called 'amr_shared' for further analysis
amr_shared <- inner_join(amr_matrix,metadata, by = c('sample' = 'renamed')) %>% 
  inner_join (metaxa_df) %>%
  inner_join (gene_length_df) %>%
  select(-c('sent_name', 'sample.y')) %>% # remove unnecessary columns
  select(c('sample','treatment','temperature','header', 'type', 'class', 'mechnism','group','hits','gene_length','num_of_16srRNA')) %>% # reorganize the columns
  mutate(normalized_arg_abundance= (hits*1432)/(gene_length*num_of_16srRNA)) # perform normalization for the arg counts per sample
# save the combined dataframe for future use
write.xlsx(amr_shared,'amr_shared.xlsx')

  
  
  
  