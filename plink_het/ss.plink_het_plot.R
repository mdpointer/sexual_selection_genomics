library(tidyverse)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)



## read in population file to add treatment info
# there is one 'wrong' name in this list, but it is also wrong in the VCF so leave it here and correct after
source("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/sample_info/generate_population.list.full.R")
population_info <- generate_ss_population.info.full.list("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/sample_info/sample_names_file_182_from_VCF.txt") %>%
  rename(ind = name) %>%
  mutate( line_rep = case_when( ind== "polym4" ~ "C",
                                TRUE ~ line_rep)) %>%
  mutate(treat_rep = paste(treatment, line_rep, sep="_"))



data_file <- "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/plink_het/ss_genomewide_heterozygosity.het"

dd_het <- read.delim(data_file, header=TRUE, sep = "") %>% 
  mutate(ind = case_when( IID == "L004" ~ paste(FID, IID, sep="_"),
                          TRUE ~ IID)) %>% 
  select(-FID, -IID) %>% 
  mutate(het_per_kb = (N.NM. - O.HOM.)/ ( 147631796 / 1000) )


dd <- full_join(dd_het, population_info) %>% 
  select(-seq_run) %>% 
  filter(treatment %in% c("GA1","mono","poly"))

F_plot <- dd %>% 
  ggplot( aes( x= treatment, y=F)) +
  geom_jitter(aes(colour=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  labs(title="Inbreeding coefficient") +
  theme(text = element_text(size = 20))




lmer


het_per_kb_plot <- dd %>% 
  ggplot( aes( x= treatment, y=het_per_kb)) +
  geom_jitter(aes(colour=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  labs(title="Hets per 1000bp") +
  theme(text = element_text(size = 20))

het_plots <- ggarrange(F_plot, het_per_kb_plot, common.legend = T)

ggsave("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/plink_het/ss.plink_heterozygosity.png",
       het_plots,
       device="png",
       units="mm",
       height=200,
       width=300)

  

dd_mod <- dd %>% 
  filter(treatment %in% c("mono","poly"))
model <- lmer(F ~ treatment + (1 | treat_rep), data = dd_mod)
summary(model)

qqnorm(resid(model)); qqline(resid(model))
plot(model)  # plots residuals and random effects

emmeans(model, pairwise ~ treatment)
confint(model)

vc <- as.data.frame(VarCorr(model))
icc <- vc$vcov[1] / sum(vc$vcov)
round(icc, 3)


dd_mod %>% 
  group_by(treatment) %>% 
  summarise(mean_F = mean(F),
            se_F = sd(F) / sqrt(n()))
