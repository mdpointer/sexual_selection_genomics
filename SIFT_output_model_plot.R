library(tidyverse)
library(ggpubr)
library(lme4)
library(lmerTest)

setwd("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/SIFT")

##################################  Get the info on all the samples
source("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/sample_info/generate_population.list.full.R")
population_info <- generate_ss_population.info.full.list("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/sample_info/sexual_selection_samplenames") %>%
  rename(sample = name) %>%
  mutate( line_rep = case_when( sample== "polym4" ~ "C",
                                TRUE ~ line_rep)) %>%
  mutate(treat_rep = paste(treatment, line_rep, sep="_"))




### BEGINNGING WITH ALL VARIANTS (IE INCL NONSENSE)

################################## ALL SNP - DELETERIOUS (MISSENSE & NONSENSE) - HIGH CONFIDENCE - AUTOSOMES

dd_all_del_highconf_autosomes <- read_table("codingRegions_all_snp_deleterious_high_conf_stats_autosomes.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="all snp",
         SIFT_impact="deleterious",
         conf="high",
         chrom="autosomes")


################################## ALL SNP - DELETERIOUS (MISSENSE & NONSENSE) - HIGH CONFIDENCE - X

dd_all_del_highconf_X <- read_table("codingRegions_all_snp_deleterious_high_conf_stats_X.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="all snp",
         SIFT_impact="deleterious",
         conf="high",
         chrom="X")




################################## ALL SNP - DELETERIOUS (MISSENSE & NONSENSE) - HIGH & LOW CONFIDENCE COMBINED - AUTOSOMES

dd_all_del_allconf_autosomes <- read_table("codingRegions_all_snp_deleterious_stats_autosomes.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="all snp",
         SIFT_impact="deleterious",
         conf="all",
         chrom="autosomes")


################################## ALL SNP - DELETERIOUS (MISSENSE & NONSENSE) - HIGH & LOW CONFIDENCE COMBINED - X

dd_all_del_allconf_X <- read_table("codingRegions_all_snp_deleterious_high_conf_stats_X.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="all snp",
         SIFT_impact="deleterious",
         conf="all",
         chrom="X")




################################## ALL SNP - TOLERATED (MISSENSE & NONSENSE) - AUTOSOMES

dd_all_tol_autosomes <- read_table("codingRegions_all_snp_tolerated_stats_autosomes.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="all snp",
         SIFT_impact="tolerated",
         conf="tol",
         chrom="autosomes")


################################## ALL SNP - TOLERATED (MISSENSE & NONSENSE) - X

dd_all_tol_X <- read_table("codingRegions_all_snp_tolerated_stats_X.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="all snp",
         SIFT_impact="tolerated",
         conf="tol",
         chrom="X")








### BEGINNGING WITH ONLY MISSENSE VARIANTS (IE EXCL NONSENSE)

################################## MISSENSE SNP - DELETERIOUS - HIGH CONFIDENCE - AUTOSOMES

dd_missense_del_highconf_autosomes <- read_table("codingRegions_missense_deleterious_high_conf_stats_autosomes.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="missense",
         SIFT_impact="deleterious",
         conf="high",
         chrom="autosomes")


################################## MISSENSE - DELETERIOUS - HIGH CONFIDENCE - X

dd_missense_del_highconf_X <- read_table("codingRegions_missense_deleterious_high_conf_snp_stats_X.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="missense",
         SIFT_impact="deleterious",
         conf="high",
         chrom="X")




################################## MISSENSE - DELETERIOUS - HIGH & LOW CONFIDENCE COMBINED - AUTOSOMES

dd_missense_del_allconf_autosomes <- read_table("codingRegions_missense_deleterious_stats_autosomes.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="missense",
         SIFT_impact="deleterious",
         conf="all",
         chrom="autosomes")


################################## MISSENSE - DELETERIOUS - HIGH & LOW CONFIDENCE COMBINED - X

dd_missense_del_allconf_X <- read_table("codingRegions_missense_deleterious_snp_stats_X.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="missense",
         SIFT_impact="deleterious",
         conf="all",
         chrom="X")


################################## MISSENSE - TOLERATED - AUTOSOMES

dd_missense_tol_autosomes <- read_table("codingRegions_missense_tolerated_stats_autosomes.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="missense",
         SIFT_impact="tolerated",
         conf="tol",
         chrom="autosomes")


################################## MISSENSE - TOLERATED - X

dd_missense_tol_X <- read_table("codingRegions_missense_tolerated_snp_stats_X.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="missense",
         SIFT_impact="tolerated",
         conf="tol",
         chrom="X")



# SYNONYMOUS SITE COUNTS
################################## SYNONYMOUS on AUTOSOMES

dd_synonymous_autosomes <- read_table("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/SIFT/codingRegions_synonymous_snp_stats_autosomes.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="synonymous",
         SIFT_impact="synonymous",
         conf="syn",
         chrom="autosomes")


################################## SYNONYMOUS on X

dd_synonymous_X <- read_table("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/SIFT/codingRegions_synonymous_snp_stats_X.txt", skip=2, col_names=F) %>% 
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
  mutate(dataset="synonymous",
         SIFT_impact="synonymous",
         conf="syn",
         chrom="X")





##################################
##################################  PULL ALL THE DATA TOGETHER
################################## 

dd <- bind_rows(
  dd_all_del_highconf_autosomes,
  dd_all_del_highconf_X,
  dd_all_del_allconf_autosomes,
  dd_all_del_allconf_X,
  dd_all_tol_autosomes,
  dd_all_tol_X,
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
  dd_missense_del_allconf_autosomes,
  dd_missense_del_allconf_X,
  dd_missense_tol_autosomes,
  dd_missense_tol_autosomes,
  dd_synonymous_autosomes,
  dd_synonymous_X
  ) %>% 
  left_join(population_info) %>%
  pivot_longer(cols=c(nRefHom, nNonRefHom, nHet), names_to="genotype", values_to="n") %>% 
  mutate(
    # SIFT_impact = case_when( SIFT_impact=="low" ~ "Low impact - tolerated",
    #                                 SIFT_impact=="moderate" ~ "Moderate impact - tolerated",
    #                                 SIFT_impact=="high" ~ "High impact - deleterious",
    #                                 SIFT_impact=="modifier" ~ "Modifier impact",
    #                                 TRUE ~ SIFT_impact),
         treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Polyandry" ))) %>% 
  filter(treatment %in% c("Ancestral", "Monandry", "Polyandry"))





################
################ PLOTS


################ PLOT: nRefHom & nHet ACROSS HIGH:LOW IMPACT & 5 SS TREATMENTS
del_vs_syn_sift_boxplot <- dd %>%
  mutate(genotype= case_when( genotype=="nRefHom" ~ "RefHom",
                              genotype=="nNonRefHom" ~ "AltHom",
                              genotype=="nHet" ~ "Het")) %>%
  mutate(genotype = factor(genotype, levels=c("RefHom", "Het", "AltHom"))) %>% 
  # mutate(SIFT_impact = factor(SIFT_impact, levels=c("Low impact - tolerated", "Moderate impact - tolerated","High impact - deleterious"))) %>% 
  filter(genotype!="RefHom") %>% 
  filter(chrom=="autosomes") %>% 
  filter(dataset %in% c("all snp", "synonymous")) %>% 
  filter(conf %in% c("high","syn")) %>%
  ggplot( aes( x= genotype, y=n, fill=genotype)) +
  geom_jitter(aes(colour=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  facet_grid(SIFT_impact~treatment, scales="free_y") +
  theme_bw() +
  labs(title="_hom.vs.het") +
  theme(text = element_text(size = 20))

ggsave("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/SIFT/del_vs_syn_sift_boxplot.png",
       del_vs_syn_sift_boxplot,
       device="png",
       units="mm",
       height=200,
       width=360)

################ PLOT Number of Homozygous deleterious loci (burden)
dd %>% 
  filter(SIFT_impact=="deleterious") %>% 
  filter(chrom=="autosomes") %>% 
  filter(genotype=="nNonRefHom") %>%
  filter(dataset=="all snp") %>% 
  filter(conf %in% c("high", "tol")) %>% 
  ggplot(aes(x=treatment, y=n)) +
  geom_jitter(aes(colour=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  # facet_wrap(~sex) +
  theme_bw() +
  labs(title="Hom deleterious 'burden")


# dd_model_nNonRefHom <- dd %>% 
#   filter(SIFT_impact=="High impact - deleterious") %>% 
#   filter(chrom=="autosomes") %>% 
#   filter(genotype=="nNonRefHom")
# 
# m1 <- lmer(n ~ relevel(treatment, ref="Female biased") + (1|treat_rep), data=dd_model_nNonRefHom)
# summary(m1)


################ PLOT Number of Homozygous deleterious loci (burden) & synonymous
dd %>% 
  filter( SIFT_impact %in% c("synonymous","deleterious")) %>%
  filter(chrom=="autosomes") %>% 
  filter(genotype != "nRefHom") %>%
  filter(dataset %in% c("all snp", "synonymous")) %>% 
  filter(conf %in% c("high", "syn")) %>% 
  ggplot(aes(x=genotype, y=n, fill=genotype)) +
  geom_jitter(aes(colour=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  facet_grid(SIFT_impact~treatment, scales="free_y") +
  theme_bw() +
  labs(title="Hom deleterious 'burden'")



################ PLOT Ratio of deleterious:synonymous sites for NonRefHom
p1 <- dd %>% 
  filter( SIFT_impact %in% c("synonymous","deleterious")) %>%
  filter(chrom=="autosomes") %>% 
  filter(genotype == "nNonRefHom") %>%
  filter(dataset %in% c("all snp", "synonymous")) %>% 
  filter(conf %in% c("high", "syn")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, n) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=n) %>% 
  mutate(ratio_del_syn=deleterious / synonymous) %>% 
  ggplot(aes(x=treatment, y=ratio_del_syn, fill=treatment)) +
  geom_jitter(aes(colour=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  # facet_wrap(~sex) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title="Ratio of deleterious:synonymous sites for NonRefHom")



################ PLOT Ratio of deleterious:synonymous sites for Het
p2 <- dd %>% 
  filter( SIFT_impact %in% c("synonymous","deleterious")) %>%
  filter(chrom=="autosomes") %>% 
  filter(genotype == "nHet") %>%
  filter(dataset %in% c("all snp", "synonymous")) %>% 
  filter(conf %in% c("high", "syn")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, n) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=n) %>% 
  mutate(ratio_del_syn=deleterious / synonymous) %>% 
  ggplot(aes(x=treatment, y=ratio_del_syn, fill=treatment)) +
  geom_jitter(aes(colour=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  # facet_wrap(~sex) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title="Ratio of deleterious:synonymous sites for Het")

ggarrange(p1,p2,common.legend = T)





## number of sites with >=1 deleterious allele
dd <- bind_rows(
  dd_all_del_highconf_autosomes,
  dd_all_del_highconf_X,
  dd_all_del_allconf_autosomes,
  dd_all_del_allconf_X,
  dd_all_tol_autosomes,
  dd_all_tol_X,
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
  dd_missense_del_allconf_autosomes,
  dd_missense_del_allconf_X,
  dd_missense_tol_autosomes,
  dd_missense_tol_autosomes,
  dd_synonymous_autosomes,
  dd_synonymous_X
) %>% 
  left_join(population_info) %>%
  filter(SIFT_impact %in% c("deleterious", "synonymous")) %>% 
  filter(chrom=="autosomes") %>% 
  filter(dataset %in% c("all snp", "synonymous")) %>% 
  filter(conf %in% c("high", "syn")) %>% 
  mutate(total_snps = nHet + nNonRefHom) %>% 
  mutate(SIFT_impact = case_when( SIFT_impact=="low" ~ "Low impact - tolerated",
                                    SIFT_impact=="moderate" ~ "Moderate impact - tolerated",
                                    SIFT_impact=="high" ~ "High impact - deleterious",
                                    SIFT_impact=="modifier" ~ "Modifier impact",
                                    TRUE ~ SIFT_impact),
         treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Enforced polyandry" ))) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=total_snps) %>% 
  mutate( ratio_del_syn_totalsnps = deleterious / synonymous) %>% 
  filter(treatment %in% c("Ancestral","Monandry","Polyandry")) %>% 
  ggplot(aes(x=treatment, y=ratio_del_syn_totalsnps, fill=treatment)) +
  geom_jitter(aes(colour=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  # facet_wrap(~sex) +
  theme_bw() +
  labs(title="Ratio of deleterious snps:syn snps present in the individual")



## number of deleterious alleles
bind_rows(
  dd_all_del_highconf_autosomes,
  dd_all_del_highconf_X,
  dd_all_del_allconf_autosomes,
  dd_all_del_allconf_X,
  dd_all_tol_autosomes,
  dd_all_tol_X,
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
  dd_missense_del_allconf_autosomes,
  dd_missense_del_allconf_X,
  dd_missense_tol_autosomes,
  dd_missense_tol_autosomes,
  dd_synonymous_autosomes,
  dd_synonymous_X
) %>% 
  left_join(population_info) %>%
  filter(SIFT_impact %in% c("deleterious", "synonymous")) %>% 
  filter(chrom=="autosomes") %>% 
  filter(dataset %in% c("all snp", "synonymous")) %>% 
  filter(conf %in% c("high", "syn")) %>% 
  mutate(total_alleles = nHet + 2*nNonRefHom) %>% 
  mutate(SIFT_impact = case_when( SIFT_impact=="low" ~ "Low impact - tolerated",
                                  SIFT_impact=="moderate" ~ "Moderate impact - tolerated",
                                  SIFT_impact=="high" ~ "High impact - deleterious",
                                  SIFT_impact=="modifier" ~ "Modifier impact",
                                  TRUE ~ SIFT_impact),
         treatment = case_when( treatment=="GA1" ~ "Outbred ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Enforced monandry",
                                treatment=="poly" ~ "Enforced polyandry"),
         treatment = factor(treatment, levels=c("Outbred ancestral","Female biased","Male biased","Enforced monandry","Enforced polyandry" ))) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_alleles) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=total_alleles) %>% 
  mutate( ratio_del_syn_totalalleles = deleterious / synonymous) %>% 
  ggplot(aes(x=treatment, y=ratio_del_syn_totalalleles, fill=treatment)) +
  geom_jitter(aes(colour=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  # facet_wrap(~sex) +
  theme_bw() +
  labs(title="Ratio of deleterious snps:syn snps present in the individual")







##################

dd <- bind_rows(
  dd_all_del_highconf_autosomes,
  dd_all_del_highconf_X,
  dd_all_del_allconf_autosomes,
  dd_all_del_allconf_X,
  dd_all_tol_autosomes,
  dd_all_tol_X,
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
  dd_missense_del_allconf_autosomes,
  dd_missense_del_allconf_X,
  dd_missense_tol_autosomes,
  dd_synonymous_autosomes,
  dd_synonymous_X
) %>% 
  left_join(population_info)  %>%  
  filter(dataset %in% c("missense","synonymous"),
         SIFT_impact %in% c("deleterious", "tolerated", "synonymous")) %>% 
  filter(chrom=="autosomes") %>%
  filter(conf %in% c("high", "tol","syn")) %>% 
  mutate(total_snps = nHet + nNonRefHom) %>%
  mutate(treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Polyandry" ))) %>% 
  filter(treatment %in% c("Ancestral", "Monandry", "Polyandry")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=total_snps) %>% 
  mutate( ratio_misdel_syn_totalsnps = deleterious / synonymous,
          ratio_mistol_syn_totalsnps = tolerated / synonymous,
          `Replicate line` = line_rep,
          Treatment = treatment)


dd %>% 
  ggplot(aes(x=treatment, y=ratio_misdel_syn_totalsnps, fill=treatment)) +
  geom_jitter(aes(shape=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("yellow2","dodgerblue4","firebrick")) +
  labs(title="Missense deleterious ",
       y=expression(N[X]~" / Syn ratio"),
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5))

dd %>% 
  ggplot(aes(x=treatment, y=ratio_mistol_syn_totalsnps, fill=treatment)) +
  geom_jitter(aes(shape=line_rep),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("yellow2","dodgerblue4","firebrick")) +
  labs(title="Missense tolerated ",
       y=expression(N[X]~" / Syn ratio"),
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5))
