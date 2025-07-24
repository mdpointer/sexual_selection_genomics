library(tidyverse)
library(ggpubr)




data_file <- "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/vcftools_diversity/genomewide_pi_summary.txt"

dd_div <- read.delim(data_file, header=TRUE, sep = "")

dd_div %>% 
  ggplot(aes(x=Treatment,y=Mean_pi)) +
  geom_jitter(aes(shape=Rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  theme_bw()


dd_div <- dd_div %>% 
  filter(Treatment %in% c("mono","poly"))
kruskal.test(Mean_pi ~ Treatment, data = dd_div)


dd_div %>% 
  group_by(Treatment) %>% 
  summarise(mean_mean_pi = mean(Mean_pi),
            se_mean_pi = sd(Mean_pi) / sqrt(n()))
