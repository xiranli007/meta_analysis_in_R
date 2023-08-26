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

# load amr matrix 
amr_shared <- read_excel('amr_shared.xlsx')
# load metadata
metadata <- read_excel('davis_water_metadata.xlsx')
# create color palatte
qual_col_pals <- brewer.pal.info[brewer.pal.info$category=="qual",]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# using vegan package to calculate beta diversity 
# the input for function avgdist needs to be a matrix with sample name as row names 
# and the other factor (e.g.sequence or amr) as column name
# create the matrix from amr_shared at group level
amr_df<-amr_shared %>%
  group_by(sample,group) %>%
  summarise(group_abund = sum(normalized_arg_abundance), .groups = 'drop') %>%
  pivot_wider(names_from = group, values_from = group_abund) %>%
  as.data.frame() 
rownames(amr_df) <- amr_df$sample #convert the sample as rowname 
amr_df <- as.matrix(amr_df[, -1]) # remove the first sample column and now we can convert to the matrix and fed into avgdist
amr_dist <- vegdist(amr_df, method ='bray')
amr_ord <- amr_dist %>% cmdscale(., eig=TRUE, add=TRUE) # perform ordination on the distance matrix
positions <- amr_ord$points  # get the x, y coordination on the pcoa plot
colnames(positions) <- c("pcoa1", "pcoa2")  # rename the coordination
percent_explained <- 100 * amr_ord$eig / sum(amr_ord$eig) # grab the x and y axis labels for the pcoa plot i.e. percentage explained
pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE) # clean up the number to make it readable as the axis titles
labels <- c(glue("PCo Axis 1 ({pretty_pe[1]}%)"),    # create the label
            glue("PCo Axis 2 ({pretty_pe[2]}%)"))

# create a function to calculate the centroids 
# the input df should have pcoa1, pcoa2 and the metadata for the sample, 'factor'is the factor used for separation in the pcoa plot
cal_centroid <- function (df, factor){
  df %>% group_by({{factor}})%>%
    summarise(pcoa1 = mean(pcoa1), pcoa2 = mean(pcoa2),.groups = 'drop')
} 

# rejoin positions with metadata
amr_beta_plot <-positions %>% as.data.frame() %>% 
  rownames_to_column('sample') %>%
  inner_join(metadata, by = c('sample' = 'renamed')) %>%
  select(-c('sample.y','sent_name'))
# separate the data and calculate the centroid points
amr_beta_plot_20 <- amr_beta_plot %>% filter(temperature == 20)
amr_20_centroid <- cal_centroid(amr_beta_plot_20, treatment)
amr_beta_plot_25 <- amr_beta_plot %>% filter(temperature == 25)
amr_25_centroid <- cal_centroid(amr_beta_plot_25, treatment)
amr_beta_plot_30 <- amr_beta_plot %>% filter(temperature == 30)
amr_30_centroid <- cal_centroid(amr_beta_plot_30, treatment)

tiff("amr_beta_20.tiff", units="in", width=6, height=6, res=600)
amr_beta_plot_20 %>% 
  ggplot(aes(x = pcoa1, y = pcoa2, fill = treatment, color = treatment))+
  geom_point()+
  labs(x=labels[1], y=labels[2])+
  theme_pubr(base_family="Times New Roman")+
  stat_ellipse(geom="polygon", aes(group = treatment),
               type = "t",level=0.8, linewidth = 0.5,alpha = 0.2)+
  geom_point(data = amr_20_centroid , size= 3, shape=21, aes(fill = treatment))+
  scale_color_manual(name = 'Treatment',values = col_vector[5:10])+
  scale_fill_manual(name = 'Treatment',values = col_vector[5:10])+
  theme (panel.background= element_rect(fill="white"),
         panel.border=element_rect(color="black",fill=NA, size=0.7),
         axis.title = element_markdown(face = 'bold'),
         legend.text = element_markdown(face = 'bold'),
         legend.title = element_markdown(face = 'bold'),
         axis.text = element_markdown(size = 10))
dev.off() 

tiff("amr_beta_25.tiff", units="in", width=6, height=6, res=600)
amr_beta_plot_25 %>% 
  ggplot(aes(x = pcoa1, y = pcoa2, fill = treatment, color = treatment))+
  geom_point()+
  labs(x=labels[1], y=labels[2])+
  theme_pubr(base_family="Times New Roman")+
  stat_ellipse(geom="polygon", aes(group = treatment),
               type = "t",level=0.8, linewidth = 0.5,alpha = 0.2)+
  geom_point(data = amr_25_centroid , size= 3, shape=21, aes(fill = treatment))+
  scale_color_manual(name = 'Treatment',values = col_vector[5:10])+
  scale_fill_manual(name = 'Treatment',values = col_vector[5:10])+
  theme (panel.background= element_rect(fill="white"),
         panel.border=element_rect(color="black",fill=NA, size=0.7),
         axis.title = element_markdown(face = 'bold'),
         legend.text = element_markdown(face = 'bold'),
         legend.title = element_markdown(face = 'bold'),
         axis.text = element_markdown(size = 10))
dev.off()

tiff("amr_beta_30.tiff", units="in", width=6, height=6, res=600)
amr_beta_plot_30 %>% 
  ggplot(aes(x = pcoa1, y = pcoa2, fill = treatment, color = treatment))+
  geom_point()+
  labs(x=labels[1], y=labels[2])+
  theme_pubr(base_family="Times New Roman")+
  stat_ellipse(geom="polygon", aes(group = treatment),
               type = "t",level=0.8, linewidth = 0.5,alpha = 0.2)+
  geom_point(data = amr_30_centroid , size= 3, shape=21, aes(fill = treatment))+
  scale_color_manual(name = 'Treatment',values = col_vector[5:10])+
  scale_fill_manual(name = 'Treatment',values = col_vector[5:10])+
  theme (panel.background= element_rect(fill="white"),
         panel.border=element_rect(color="black",fill=NA, size=0.7),
         axis.title = element_markdown(face = 'bold'),
         legend.text = element_markdown(face = 'bold'),
         legend.title = element_markdown(face = 'bold'),
         axis.text = element_markdown(size = 10))
dev.off() 

# perform statistical analysis
# write a function to extract desired data from the distance matrix
# import the distance for all samples as dist, use f as factor to select a subset of sample
select_dist <- function (dist,f){
  sub_dist <- dist %>% as.matrix() %>% 
    as.data.frame() %>% 
    rownames_to_column('sample') %>%
    inner_join(metadata, by=c('sample'='renamed')) %>% 
    filter(temperature=={{f}}) %>%
    select (-c(sent_name, sample.y, treatment,temperature)) %>% 
    column_to_rownames('sample')
  sub_dist <- sub_dist[,colnames(sub_dist) %in% rownames(sub_dist)] %>% as.dist()
}
dist_20 <-select_dist(amr_dist,20)
dist_25 <-select_dist(amr_dist,25)
dist_30 <-select_dist(amr_dist,30)

# use adonis and pairwise adonis from FSA package to perform statistical analysis
set.seed(199701) # set a seed to obtain the same results for the permutation test
permutations = 1000

beta_20_adonis <- adonis2(dist_20 ~ treatment, amr_beta_plot_20, permutations = permutations) # p=0.001998002
beta_20_posthoc <- pairwise.adonis(dist_20,amr_beta_plot_20$treatment, p.adjust.m = 'BH') # 0.033 for each pair

beta_25_adonis <- adonis2(dist_25 ~ treatment, amr_beta_plot_25, permutations = permutations) # p=0.03996004
beta_25_posthoc <- pairwise.adonis(dist_25,amr_beta_plot_25$treatment, p.adjust.m = 'BH') # no

beta_30_adonis <- adonis2(dist_30 ~ treatment, amr_beta_plot_30, permutations = permutations) # p=0.002997003
beta_30_posthoc <- pairwise.adonis(dist_30,amr_beta_plot_30$treatment, p.adjust.m = 'BH') # no

  







