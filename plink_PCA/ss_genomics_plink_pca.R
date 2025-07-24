library(tidyverse)
library(stringr)
library(patchwork)


## read in population file to add treatment info
# there is one 'wrong' name in this list, but it is also wrong in the VCF so leave it here and correct after
source("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/sample_info/generate_population.list.full.R")
population_info <- generate_ss_population.info.full.list("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/sample_info/sample_names_file_182_from_VCF.txt") %>%
  rename(ind = name) %>%
  mutate( line_rep = case_when( ind== "polym4" ~ "C",
          TRUE ~ line_rep)) %>%
  mutate(treat_rep = paste(treatment, line_rep, sep="_")) %>% 
  filter(treatment %in% c("mono","poly", "GA1"))


# read in data
pca <- read_table("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/plink_pca/108samples/ss_snps_182samples_plink_MAF_ldpruned_50_10_0.1_PCA.eigenvec", col_names = FALSE)
eigenval <- scan("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/plink_pca/108samples/ss_snps_182samples_plink_MAF_ldpruned_50_10_0.1_PCA.eigenval")


# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
# set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))



## add info from population info to pca data
# need to add treatment, sex, batch at least
pca <- full_join(pca, population_info) %>%
  as_tibble() %>% 
  drop_na(treatment, sex) %>% 
  mutate(Sex=sex)

## create a tibble that takes each eigenval and converts to percentage variance explained
eigenval <- tibble(PC = 1:20, explained = eigenval/sum(eigenval)*100)
plot_eigenval <- eigenval %>%
  ggplot(aes(x=PC, y=explained)) +
  geom_bar(stat="identity") +
  theme_bw() +
  labs(y="% var explained")



# PLOT THE ACTUAL PCA PLOT
pca %>%
  ggplot( aes( x=PC1, y=PC2, col=treatment, shape=sex)) +
  geom_point(size=3) +
  theme_bw()

pca %>%
  ggplot( aes( x=PC3, y=PC4, col=treatment, shape=as.factor(sex))) +
  geom_point(size=3) +
  theme_bw()

pca %>%
  ggplot( aes( x=PC2, y=PC4, col=treatment, shape=sex)) +
  geom_point() +
  theme_bw()

pca %>%
  ggplot( aes( x=PC2, y=PC5, col=treatment, shape=sex)) +
  geom_point() +
  theme_bw()


pca %>%
  ggplot( aes( x=PC2, y=PC6, col=treatment, shape=sex)) +
  geom_point() +
  theme_bw()

pca %>%
  ggplot( aes( x=PC2, y=PC7, col=treatment, shape=sex)) +
  geom_point() +
  theme_bw()

pca %>%
  ggplot( aes( x=PC2, y=PC8, col=treatment, shape=sex)) +
  geom_point() +
  theme_bw()

pca %>% filter(treatment == "mono") %>%
  ggplot( aes( x=PC2, y=PC3, col=treat_rep, shape=treatment)) +
  geom_point() +
  theme_bw()




pca %>%
  ggplot( aes( x=PC1, y=PC2, col=treatment, shape=line_rep)) +
  geom_point(size=3) +
  theme_bw()

pca %>%
  ggplot( aes( x=PC3, y=PC4, col=treatment, shape=line_rep)) +
  geom_point(size=3) +
  theme_bw()








################
################


# PLOT THE ACTUAL PCA PLOT
PC1_v_2 <- pca %>%
  rename(`SS regime` = treatment) %>% 
  rename(Replicate = line_rep) %>% 
  mutate(`SS regime` = case_when(`SS regime`=="mono"~"Monandry",
                                 `SS regime`=="poly"~"Polyandry",
                                 `SS regime`=="GA1"~"Ancestral"),
         `SS regime` = factor(`SS regime`, levels = c(
           "Ancestral", "Monandry", "Polyandry" ))) %>% 
  ggplot( aes( x=PC1, y=PC2, col=`SS regime`, shape=Replicate)) + # fill= as.factor(rep)
  geom_point(size=3, alpha=0.7) +  #, shape=21, stroke =2 
  theme_bw() +
  scale_colour_manual(values=c("yellow3","dodgerblue4","firebrick")) +
  theme(text = element_text(size=13)) +
  labs(x="PC1 - 20.1%",
       y="PC2 - 14.2%")

PC3_v_4 <- pca %>% 
  rename(`SS regime` = treatment) %>% 
  rename(Replicate = line_rep) %>% 
  mutate(`SS regime` = case_when(`SS regime`=="mono"~"Monandry",
                                 `SS regime`=="poly"~"Polyandry",
                                 `SS regime`=="GA1"~"Ancestral"),
         `SS regime` = factor(`SS regime`, levels = c(
           "Ancestral", "Monandry", "Polyandry" ))) %>% 
  ggplot( aes( x=PC3, y=PC4, col=`SS regime`, shape=Replicate)) + # fill= as.factor(rep)
  geom_point(size=3, alpha=0.7) + # , shape=21, stroke =2
  theme_bw() +
  scale_colour_manual(values=c("yellow3","dodgerblue4","firebrick")) +
  theme(text = element_text(size=13)) +
  labs(x="PC3 - 11.8%",
       y="PC4 - 9.2%")


PC1_v_2_sex <- pca %>%
  rename(`SS regime` = treatment) %>% 
  mutate(`SS regime` = case_when(`SS regime`=="mono"~"Monandry",
                                 `SS regime`=="poly"~"Polyandry",
                                 `SS regime`=="GA1"~"Ancestral"),
         `SS regime` = factor(`SS regime`, levels = c(
           "Ancestral", "Monandry", "Polyandry" ))) %>% 
  ggplot( aes( x=PC1, y=PC2, col=`SS regime`, shape=Sex)) +
  geom_point(size=3, alpha=0.7) +
  theme_bw() +
  scale_colour_manual(values=c("yellow3","dodgerblue4","firebrick")) +
  theme(text = element_text(size=13)) +
  labs(x="PC1 - 20.1%",
       y="PC2 - 14.2%")


legend <- get_legend(PC1_v_2, position="bottom")

out_plot_1 <- ggarrange(
  PC1_v_2,
  PC3_v_4,
  #PC1_v_2_sex,
  common.legend = T,
  labels="AUTO",
  # legend.grob = legend,
  legend="bottom"
)

# out_plot_2 <- ggarrange(
#   plot_eigenval,
#   PC1_v_2_sex,
#   common.legend = T,
#   labels="AUTO",
#   legend="bottom"
# )



ggsave(
  filename = "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/plink_pca/disp_PCA_2.png",
  plot = PC1_v_2_sex,
  device = "png",
  units="mm",
  height=100,
  width=120,
  dpi=300
)



