library(tidyverse)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)

setwd("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/snpEff")

##################################  Get the info on all the samples
source("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/sample_info/generate_population.list.full.R")
population_info <- generate_ss_population.info.full.list("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/sample_info/sexual_selection_samplenames") %>%
  rename(sample = name) %>%
  mutate( line_rep = case_when( sample== "polym4" ~ "C",
                                TRUE ~ line_rep)) %>%
  mutate(treat_rep = paste(treatment, line_rep, sep="_"))




################################## HIGH IMPACT on AUTOSOMES

dd_h_a <- read_table("codingRegions_high_impact_stats_autosomes.txt", skip=2, col_names=F) %>% 
  select(3:6) %>% 
  rename(sample=X3,
         nRefHom=X4,
         nNonRefHom=X5,
         nHet=X6
         ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="high",
         chrom="autosomes")


################################## HIGH IMPACT on X

dd_h_x <- read_table("codingRegions_high_impact_stats_X.txt", skip=2, col_names=F) %>% 
  rename(sample=X1,
         nRefHom=X2,
         nNonRefHom=X3,
         nHet=X4
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="high",
         chrom="X")



################################## MODERATE IMPACT on AUTOSOMES

dd_m_a <- read_table("codingRegions_moderate_impact_stats_autosomes.txt", skip=2, col_names=F) %>% 
  select(3:6) %>% 
  rename(sample=X3,
         nRefHom=X4,
         nNonRefHom=X5,
         nHet=X6
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="moderate",
         chrom="autosomes")



################################## MODERATE IMPACT on X

dd_m_x <- read_table("codingRegions_moderate_impact_stats_X.txt", skip=2, col_names=F) %>% 
  rename(sample=X1,
         nRefHom=X2,
         nNonRefHom=X3,
         nHet=X4
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="moderate",
         chrom="X")




################################## LOW IMPACT on AUTOSOMES

dd_l_a <- read_table("codingRegions_low_impact_stats_autosomes.txt", skip=2, col_names=F) %>% 
  select(3:6) %>% 
  rename(sample=X3,
         nRefHom=X4,
         nNonRefHom=X5,
         nHet=X6
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="low",
         chrom="autosomes")




################################## LOW IMPACT on X

dd_l_x <- read_table("codingRegions_low_impact_stats_X.txt", skip=2, col_names=F) %>% 
  rename(sample=X1,
         nRefHom=X2,
         nNonRefHom=X3,
         nHet=X4
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="low",
         chrom="X")



################################## MODIFIER IMPACT on AUTOSOMES

dd_modifier_a <- read_table("codingRegions_modifier_impact_stats_autosomes.txt", skip=2, col_names=F) %>% 
  select(3:6) %>% 
  rename(sample=X3,
         nRefHom=X4,
         nNonRefHom=X5,
         nHet=X6
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="modifier",
         chrom="autosomes")




################################## MODIFIER IMPACT on X

dd_modifier_x <- read_table("codingRegions_high_impact_stats_X.txt", skip=2, col_names=F) %>% 
  rename(sample=X1,
         nRefHom=X2,
         nNonRefHom=X3,
         nHet=X4
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="modifier",
         chrom="X")



################################## ALL SNP on AUTOSOMES

dd_allsnp_autosomes <- read_table("codingRegions_all_snp_stats_autosomes.txt", skip=2, col_names=F) %>% 
  select(3:6) %>% 
  rename(sample=X3,
         nRefHom=X4,
         nNonRefHom=X5,
         nHet=X6
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="allsnp",
         chrom="autosomes")


################################## ALL SNP on X

dd_allsnp_x <- read_table("codingRegions_all_snp_stats_x.txt", skip=2, col_names=F) %>% 
  rename(sample=X1,
         nRefHom=X2,
         nNonRefHom=X3,
         nHet=X4
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="allsnp",
         chrom="X")

################################## NONSENSE on AUTOSOMES

dd_nonsense_autosomes <- read_table("codingRegions_nonsense_snp_stats_autosomes.txt", skip=2, col_names=F) %>% 
  select(3:6) %>% 
  rename(sample=X3,
         nRefHom=X4,
         nNonRefHom=X5,
         nHet=X6
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="nonsense",
         chrom="autosomes")


################################## NONSENSE on X

dd_nonsense_x <- read_table("codingRegions_nonsense_snp_stats_x.txt", skip=2, col_names=F) %>% 
  rename(sample=X1,
         nRefHom=X2,
         nNonRefHom=X3,
         nHet=X4
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="nonsense",
         chrom="X")


################################## SYNONYMOUS on AUTOSOMES

dd_synonymous_autosomes <- read_table("codingRegions_synonymous_snp_stats_autosomes.txt", skip=2, col_names=F) %>% 
  select(3:6) %>% 
  rename(sample=X3,
         nRefHom=X4,
         nNonRefHom=X5,
         nHet=X6
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="synonymous",
         chrom="autosomes")


################################## SYNONYMOUS on X

dd_synonymous_x <- read_table("codingRegions_synonymous_snp_stats_x.txt", skip=2, col_names=F) %>% 
  rename(sample=X1,
         nRefHom=X2,
         nNonRefHom=X3,
         nHet=X4
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="synonymous",
         chrom="X")


################################## MISSENSE on AUTOSOMES

dd_missense_autosomes <- read_table("codingRegions_missense_snp_stats_autosomes.txt", skip=2, col_names=F) %>% 
  select(3:6) %>% 
  rename(sample=X3,
         nRefHom=X4,
         nNonRefHom=X5,
         nHet=X6
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="missense",
         chrom="autosomes")

################################## MISSENSE on X

dd_missense_x <- read_table("codingRegions_missense_snp_stats_x.txt", skip=2, col_names=F) %>% 
  rename(sample=X1,
         nRefHom=X2,
         nNonRefHom=X3,
         nHet=X4
  ) %>% 
  mutate(sample = case_when(
    grepl("_L004$", sample) ~ sub("_L004$", "", sample),
    TRUE ~ sample
  ),
  sample = case_when( sample=="polym4" ~"polyCm4",
                      TRUE ~ sample)) %>% 
  mutate(snpEff_impact="missense",
         chrom="X")

##################################
##################################  PULL ALL THE DATA TOGETHER
################################## 

dd <- bind_rows(dd_h_a, dd_h_x, dd_m_a, dd_m_x, dd_l_a, dd_l_x, dd_modifier_a, dd_modifier_x,
                dd_allsnp_autosomes, dd_allsnp_x, dd_nonsense_autosomes, dd_nonsense_x, dd_synonymous_autosomes,
                dd_synonymous_x, dd_missense_autosomes, dd_missense_x) %>% 
  left_join(population_info) %>%
  pivot_longer(cols=c(nRefHom, nNonRefHom, nHet), names_to="genotype", values_to="n") %>% 
  mutate(snpEff_impact = case_when( snpEff_impact=="low" ~ "Low impact",
                                    snpEff_impact=="moderate" ~ "Moderate impact",
                                    snpEff_impact=="high" ~ "High impact",
                                    snpEff_impact=="modifier" ~ "Modifier",
                                    TRUE ~ snpEff_impact),
         treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Polyandry" ))) %>% 
  filter(treatment %in% c("Ancestral", "Monandry", "Polyandry"))



################
################ PLOTS


################ PLOT Number of Homozygous deleterious loci (burden) & synonymous
A <- dd %>% 
  filter( snpEff_impact %in% c("synonymous","High impact")) %>%
  filter(chrom=="autosomes") %>% 
  filter(genotype != "nRefHom") %>%
  mutate(genotype=case_when(genotype=="nHet" ~ "Heterzygote",
                            genotype=="nNonRefHom" ~ "Homozygote"),
         `Replicate line` = line_rep) %>% 
  ggplot(aes(x=genotype, y=n, fill=genotype)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  facet_grid(snpEff_impact~treatment, scales="free_y") +
  theme_bw() +
  labs(title="A",
       y="Number of variants",
       x="Genotype") +
  scale_fill_manual(values=c("darkgoldenrod","darkgreen")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 10))

################ PLOT Number of Homozygous deleterious loci (burden)
B <- dd %>%
  filter(snpEff_impact=="High impact") %>% 
  filter(chrom=="autosomes") %>% 
  filter(genotype=="nNonRefHom") %>%
  ggplot(aes(x=treatment, y=n)) +
  geom_jitter(aes(shape=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5, fill="darkgreen")+
  # facet_wrap(~sex) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 10)) +
  labs(title="B",
       y="Hom deleterious variants per indiv",
       x="Sexual selection regime")

out_plot_snpEff <- ggarrange(A,B, common.legend = T, legend = "bottom")





dd_model_nNonRefHom <- dd %>%
  filter(snpEff_impact=="High impact") %>%
  filter(chrom=="autosomes") %>%
  filter(genotype=="nNonRefHom")

m_high_impact <- lmer(n ~ treatment + (1|treat_rep), data=dd_model_nNonRefHom)
summary(m_high_impact)

# Compute estimated marginal means (EMMs)
emm <- emmeans(m_high_impact, pairwise ~ treatment)

# Pairwise comparisons between treatments
pairwise_results <- pairs(emm)
summary(pairwise_results)


################ PLOT Ratio of deleterious:synonymous sites for NonRefHom
p1 <- dd %>% 
  filter( snpEff_impact %in% c("synonymous","High impact - deleterious")) %>%
  filter(chrom=="autosomes") %>% 
  filter(genotype == "nNonRefHom") %>%
  pivot_wider(names_from=snpEff_impact, values_from=n) %>% 
  mutate(ratio_del_syn=`High impact - deleterious` / synonymous) %>% 
  ggplot(aes(x=treatment, y=ratio_del_syn, fill=treatment)) +
  geom_jitter(aes(colour=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  # facet_wrap(~sex) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title="Ratio of deleterious:synonymous sites for NonRefHom")



################ PLOT Ratio of deleterious:synonymous sites for Het
p2 <- dd %>% 
  filter( snpEff_impact %in% c("synonymous","High impact - deleterious")) %>%
  filter(chrom=="autosomes") %>% 
  filter(genotype == "nHet") %>%
  pivot_wider(names_from=snpEff_impact, values_from=n) %>% 
  mutate(ratio_del_syn=`High impact - deleterious` / synonymous) %>% 
  ggplot(aes(x=treatment, y=ratio_del_syn, fill=treatment)) +
  geom_jitter(aes(colour=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  # facet_wrap(~sex) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title="Ratio of deleterious:synonymous sites for Het")

ggarrange(p1,p2,common.legend = T)





dd <- bind_rows(dd_h_a, dd_h_x, dd_m_a, dd_m_x, dd_l_a, dd_l_x, dd_modifier_a, dd_modifier_x,
          dd_allsnp_autosomes, dd_allsnp_x, dd_nonsense_autosomes, dd_nonsense_x, dd_synonymous_autosomes,
          dd_synonymous_x, dd_missense_autosomes, dd_missense_x) %>% 
  left_join(population_info) %>% 
  filter(snpEff_impact %in% c("high", "synonymous", "moderate", "low")) %>% 
  filter(chrom=="autosomes") %>% 
  mutate(total_snps = nHet + nNonRefHom) %>% 
  mutate(snpEff_impact = case_when( snpEff_impact=="low" ~ "Low impact",
                                    snpEff_impact=="moderate" ~ "Moderate impact",
                                    snpEff_impact=="high" ~ "High impact",
                                    TRUE ~ snpEff_impact),
         treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Polyandry" ))) %>% 
  filter(treatment %in% c("Ancestral", "Monandry", "Polyandry")) %>% 
  select(sample, snpEff_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=snpEff_impact, values_from=total_snps) %>% 
  mutate( ratio_del_syn_totalsnps = `High impact` / synonymous,  ## !!!
          ratio_mod_syn_totalsnps = `Moderate impact` / synonymous,
          ratio_low_syn_totalsnps = `Low impact` / synonymous,
          `Replicate line` = line_rep,
          Treatment = treatment)



A <- dd %>% 
    ggplot(aes(x=Treatment, y=ratio_del_syn_totalsnps, fill=Treatment)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("yellow2","dodgerblue4","firebrick")) +
  labs(title="High impact",
       y=expression(N[X]~" / Syn ratio"),
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5))
  
B <- dd %>% 
  ggplot(aes(x=Treatment, y=ratio_mod_syn_totalsnps, fill=Treatment)) +
  geom_jitter(aes(shape=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("yellow2","dodgerblue4","firebrick")) +
  labs(title="Moderate impact",
       y=expression(N[X]~" / Syn ratio"),
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5))



out_plot_snpEff_impact_syn_ratio <- ggarrange(A,B,common.legend = T, legend = "bottom", labels = "AUTO")



out_plot_snpEff_impact_syn_ratio_sex <- dd %>% 
  ggplot(aes(x=sex, y=ratio_del_syn_totalsnps, fill=sex)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  # facet_wrap(~treatment)+
  scale_fill_manual(values=c("lightpink3","slategray3")) +
  labs(title="High impact",
       y=expression(N[X]~" / Syn ratio"),
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/snpEff/nSNPs_impact_syn_ratio_sex.png",
       out_plot_snpEff_impact_syn_ratio_sex,
       device="png",
       units="mm",
       width=120,
       height=100)





####################

m_high_impact_int <- lmer(ratio_del_syn_totalsnps ~ treatment + sex + treatment*sex +(1|line_rep), data=dd)
summary(m_high_impact_int)

m_high_impact <- lmer(ratio_del_syn_totalsnps ~ treatment + sex + treatment*sex +(1|line_rep), data=dd)
summary(m_high_impact)

# Compute estimated marginal means (EMMs)
emm_high <- emmeans(m_high_impact, pairwise ~ treatment)
emm_high
####################

m_mod_impact_int <- lmer(ratio_mod_syn_totalsnps ~ treatment + sex + treatment*sex + (1|line_rep), data=dd)
summary(m_mod_impact_int)

m_mod_impact <- lmer(ratio_mod_syn_totalsnps ~ treatment + sex + treatment*sex + (1|line_rep), data=dd)
summary(m_mod_impact)

# Compute estimated marginal means (EMMs)
emm_mod <- emmeans(m_mod_impact, pairwise ~ treatment)
emm_mod

####################



dd %>%
  ggplot(aes(ratio_del_syn_totalsnps)) +
  geom_histogram() +
  theme_bw()

dd %>%
  ggplot(aes(ratio_mod_syn_totalsnps)) +
  geom_histogram() +
  theme_bw()
  











dd <- bind_rows(dd_h_a, dd_h_x, dd_m_a, dd_m_x, dd_l_a, dd_l_x, dd_modifier_a, dd_modifier_x,
                dd_allsnp_autosomes, dd_allsnp_x, dd_nonsense_autosomes, dd_nonsense_x, dd_synonymous_autosomes,
                dd_synonymous_x, dd_missense_autosomes, dd_missense_x) %>% 
  left_join(population_info) %>% 
  filter(snpEff_impact %in% c("nonsense", "synonymous")) %>% 
  filter(chrom=="autosomes") %>% 
  mutate(total_snps = nHet + nNonRefHom) %>% 
  mutate(snpEff_impact = case_when( snpEff_impact=="nonsense" ~ "nonsense",
                                    TRUE ~ snpEff_impact),
         treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Polyandry" ))) %>% 
  filter(treatment %in% c("Ancestral", "Monandry", "Polyandry")) %>% 
  select(sample, snpEff_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=snpEff_impact, values_from=total_snps) %>% 
  mutate( ratio_nonsense_syn_totalsnps = nonsense / synonymous,
          `Replicate line` = line_rep,
          Treatment = treatment)


dd %>% 
  ggplot(aes(x=Treatment, y=ratio_nonsense_syn_totalsnps, fill=Treatment)) +
  geom_jitter(aes(shape=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("yellow2","dodgerblue4","firebrick")) +
  labs(title="Moderate impact",
       y=expression(N[X]~" / Syn ratio"),
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5))
