library(tidyverse)
library(readxl)
library(openxlsx)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(ggtext)
library (glue)
library(ggpubr)
library(ggh4x)

# load amr_shared df generate from AMR_analysis.R
amr_shared <- read_excel('amr_shared.xlsx')
# make composition plot on amr class level
# create a summary table for samples
amr_plot <- amr_shared %>%
  mutate(class = case_when(class == 'Biocide_and_metal_resistance' ~ 'Biocide and metal resistance',
                           class == 'Mercury_resistance' ~ 'Metal resistance',
                           class == 'Multi-biocide_resistance' ~ 'Multi-biocide resistance',
                           class == 'Multi-drug_resistance' ~ 'Multi-drug resistance',
                           class == 'Multi-metal_resistance' ~ 'Multi-metal resistance',
                           class == 'Peroxide_resistance' ~ 'Biocide resistance',
                           class == 'Chromium_resistance' ~'Metal resistance',
                           class == 'Copper_resistance' ~ 'Metal resistance',
                           class == 'Cationic_antimicrobial_peptides' ~ 'CAPs',
                           class == 'Drug_and_biocide_resistance' ~ 'Drug and biocide resistance',
                           class == 'Mycobacterium_tuberculosis-specific_Drug' ~ 'MTSD', 
                           TRUE ~ class)) %>% 
  group_by(sample) %>% 
  mutate(total_arg = sum(normalized_arg_abundance)) %>%
  ungroup() %>%
  group_by(sample,group) %>%
  mutate (gene_rel_abund = 100* (sum(normalized_arg_abundance)/total_arg)) %>%
  ungroup() %>% 
  group_by (sample, class) %>%
  mutate(class_rel_abund = 100*(sum(normalized_arg_abundance)/total_arg)) %>% 
  ungroup()

amr_gene <- distinct(amr_plot, sample, group, .keep_all= TRUE) %>%
  select(c('sample','treatment','temperature','group','gene_rel_abund'))
amr_class <- distinct (amr_plot, sample, class, .keep_all = TRUE) %>%
  select(c('sample','treatment','temperature','class','class_rel_abund'))%>%
  mutate(num = strtoi(str_extract(sample,"\\d+")))
#write.xlsx(amr_class,'amr_class.xlsx')
#write.xlsx(amr_gene,'amr_class.xlsx')
# create color palette
qual_col_pals <- brewer.pal.info[brewer.pal.info$category=="qual",]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# make composition plot on class level
tiff("amr_class.tiff", units="in", width=12, height=10, res=600)
amr_class %>%
  ggplot(aes(x = sample, y = class_rel_abund, fill = class))+
  geom_col()+
  scale_y_continuous(expand=c(0,0))+
  theme_pubr()+
  facet_nested(cols = vars(temperature,treatment), scales = "free_x", space ="free")+
  scale_fill_manual(name = 'AMR class', breaks = waiver(), values = col_vector[-2])+
  labs(x='',y='Relative abundance')+
  theme(text = element_text(family = "Times New Roman", face = 'bold', color = 'black'),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_markdown(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.spacing.x = unit(-0.1,"lines"),
        panel.background= element_rect(fill="white"),
        panel.border=element_rect(color="black",fill=NA, size=1),
        strip.text = element_text(family = "Times New Roman", face="bold",size=12))
dev.off()

# make composition plot on gene level
amr_gene_pool <- amr_gene %>%
  group_by(group) %>%
  summarise(pool = max(gene_rel_abund) < 5, .groups = 'drop')
my_order <- c('floR','A16S','CAP16S','copA','GYRA','MLS23S','rpoB','rpoC','rpsA','rpsL','rrsC','rrsH','sodB','sul1','tet16S','tetG','tetR','TUFAB','Others')
tiff("amr_gene.tiff", units="in", width=12, height=10, res=600)
amr_gene %>%
  inner_join(amr_gene_pool,by = 'group') %>%
  mutate (group = if_else(pool,'Others',group)) %>%
  mutate (group = case_when(group == 'FLOR' ~ 'floR', # florfenicol resistance
                            group == 'COPA' ~ 'copA', # copper resistance
                            group == 'TETR' ~ 'tetR', # tetracycline resistance
                            group == 'TETG' ~ 'tetG',
                            group == 'RPOB' ~ 'rpoB', # rif resistsance
                            group == 'RPSA' ~ 'rpsA', # pyrazinamide(antituberculosis) resistance)
                            group == 'RPSL' ~ 'rpsL', # streptomycin resistance
                            group == 'RRSC' ~ 'rrsC', # kasugamicin resistance in Ecoli
                            group == 'RRSH' ~ 'rrsH', # spectinomycin(aminoglycoside) resistane in ecoli
                            group == 'SODB' ~ 'sodB', # oxidactive stress (superoxide dismutase) resistance
                            group == 'TET16S' ~ 'tet16S',
                            group == 'RPOCR' ~ 'rpoC', #  rif resistsance
                            group == 'SULI' ~ 'sul1',# plasmid-borne sulfonamide resistance
                            TRUE ~ group)) %>%
  mutate(group = factor(group, levels = my_order)) %>%
  #factor(group, levels = my_order) %>%
  ggplot(aes(x = sample, y = gene_rel_abund, fill = group))+
  geom_col()+
  scale_y_continuous(expand=c(0,0))+
  theme_pubr()+
  facet_nested(cols = vars(temperature,treatment), scales = "free_x", space ="free")+
  scale_fill_manual(name = 'AMR genes', breaks = waiver(), values = c(col_vector[1:18],'grey'))+
  labs(x='',y='Relative abundance')+
  theme(text = element_text(family = "Times New Roman", face = 'bold', color = 'black'),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_markdown(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.spacing.x = unit(-0.1,"lines"),
        panel.background= element_rect(fill="white"),
        panel.border=element_rect(color="black",fill=NA, size=1),
        strip.text = element_text(family = "Times New Roman", face="bold",size=12))
dev.off()








