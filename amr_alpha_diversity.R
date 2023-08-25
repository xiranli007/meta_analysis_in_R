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
amr_alpha <- amr_shared %>%
  group_by(sample) %>%
  summarise (sobs = specnumber(normalized_arg_abundance),
            shannon = diversity(normalized_arg_abundance, index = 'shannon'))%>% 
  inner_join(metadata, by = c('sample' = 'renamed')) %>%
  select(-c('sample.y', 'sent_name'))
#calculate p values
kruskal.test(sobs ~ treatment,amr_alpha, subset = temperature == 20)
kruskal.test(sobs ~ treatment,amr_alpha, subset = temperature == 25)
kruskal.test(sobs ~ treatment,amr_alpha, subset = temperature == 30)
# create labels for p values
p_values_sobs <- data.frame(
  label = c('P = 0.292','P = 0.025', 'P = 0.023'),
  temperature = c(20,25,30),
  treatment = c('CON','FLO','FLO_EW')
)
tiff("amr_sobs.tiff", units="in", width=6, height=6, res=600)
amr_alpha %>%
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
        strip.text = element_text(family = "Times New Roman", face="bold",size=12))+
  geom_text(
    data = p_values_sobs,
    mapping = (aes(x = Inf, y = -Inf, label = label)),
    hjust = 1.1,
    vjust = -0.3,
    size = 4, 
    family = 'Times New Roman',
    fontface = 2 # 2(bold), 3(italic), 4(bold.italic)
  )
dev.off()

# posthoc analysis if necessary
sobs_25 <- amr_alpha %>% filter(temperature == 25)
dunnTest(sobs_25$sobs, sobs_25$treatment,method="bh")$res %>% filter(P.adj<0.05) # CON - FLO 0.028
sobs_30 <- amr_alpha %>% filter(temperature == 30)
dunnTest(sobs_30$sobs, sobs_30$treatment,method="bh")$res %>% filter(P.adj<0.05) # CON - FLO_EW 0.018

# Now do exactly the same thing for shannon diveristy
kruskal.test(shannon ~ treatment,amr_alpha, subset = temperature == 20)
kruskal.test(shannon ~ treatment,amr_alpha, subset = temperature == 25)
kruskal.test(shannon ~ treatment,amr_alpha, subset = temperature == 30)
# create labels for p values
p_values_shannon <- data.frame(
  label = c('P = 0.211','P = 0.035', 'P =  0.092'),
  temperature = c(20,25,30),
  treatment = c('CON','FLO','FLO_EW')
)
tiff("amr_shannon.tiff", units="in", width=6, height=6, res=600)
amr_alpha %>%
  ggplot(aes(x = treatment,y = shannon, fill = treatment))+
  geom_boxplot(outlier.shape = NA, size = 0.5, width = 0.5)+
  facet_grid(~temperature, scales = "free_x", space ="free")+
  #geom_jitter(width=0.25,shape=21, color='black')+
  theme_pubr()+
  labs(x= '',y ='Shannon diversity')+
  scale_fill_manual(name = 'Treatment', values = col_vector)+
  theme(text = element_text(family = "Times New Roman", face = 'bold', color = 'black'),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.spacing = unit(-0.1,'lines'),
        panel.border = element_rect(color="black",fill=NA, size=0.7),
        strip.text = element_text(family = "Times New Roman", face="bold",size=12))+
  geom_text(
    data = p_values_shannon,
    mapping = (aes(x = Inf, y = -Inf, label = label)),
    hjust = 1.1,
    vjust = -0.3,
    size = 4, 
    family = 'Times New Roman',
    fontface = 2 # 2(bold), 3(italic), 4(bold.italic)
  )
dev.off()

# posthoc analysis if necessary
shannon_25 <- amr_alpha %>% filter(temperature == 25)
dunnTest(shannon_25$shannon, shannon_25$treatment,method="bh")$res %>% filter(P.adj<0.05) # CON - FLO_EW p.adj=0.043







