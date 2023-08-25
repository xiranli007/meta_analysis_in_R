library(tidyverse)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(ggtext)
library(stringr)
library(vegan)
library(FSA)
library (glue)
library(pairwiseAdonis)
library(ggpubr)

set.seed(199701)
# laod dataset
amr_shared <- read_excel('amr_shared.xlsx')
# load metadata
metadata <- read_excel('davis_water_metadata.xlsx')
#create color palatte
qual_col_pals <- brewer.pal.info[brewer.pal.info$category=="qual",]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# alpha diversity box plot
tiff()
amr_shared %>%
  group_by(sample) %>%
  summarise (sobs = specnumber(normalized_arg_abundance),
            shannon = diversity(normalized_arg_abundance, index = 'shannon'))%>% 
  inner_join(metadata, by = c('sample' = 'renamed')) %>%
  select(-c('sample.y', 'sent_name')) %>%
  ggplot(aes(x = treatment,y = sobs, fill = treatment))+
  geom_boxplot(outlier.shape = NA, size = 0.5, width = 0.5)+
  facet_grid(~temperature, scales = "free_x", space ="free")+
  #geom_jitter(width=0.25,shape=21, color='black')+
  theme_pubr()+
  labs(x= '',y ='Observed richness')+
  scale_fill_manual(name = 'Treatment', values = col_vector)+
  theme(text = element_text(family = "Times New Roman", face = 'bold', color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(-0.1,'lines'),
        panel.border = element_rect(color="black",fill=NA, size=0.7),
        strip.text = element_text(family = "Times New Roman", face="bold",size=12))
    










