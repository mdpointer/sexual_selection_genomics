library(tidyverse)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggsignif)
library(rlang)


#######################     GET THE RELEVANT SNPEFF DATA & PLOT - NONSENSE
#######################

setwd("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/snpEff")

##################################  Get the info on all the samples
source("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/sample_info/generate_population.list.full.R")
population_info <- generate_ss_population.info.full.list("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/sample_info/sexual_selection_samplenames") %>%
  rename(sample = name) %>%
  mutate( line_rep = case_when( sample== "polym4" ~ "C",
                                TRUE ~ line_rep)) %>%
  mutate(treat_rep = paste(treatment, line_rep, sep="_"))



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



dd_snpeff_autosomes <- bind_rows(dd_nonsense_autosomes, dd_nonsense_x, dd_synonymous_autosomes,
                                 dd_synonymous_x) %>% 
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
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(-chrom, -nRefHom, -seq_run, -treat_rep) %>% 
  pivot_wider(names_from=snpEff_impact, values_from=c(total_snps, nHet, nNonRefHom)) %>% 
  mutate( ratio_nonsense_syn_totalsnps = total_snps_nonsense / total_snps_synonymous,
          ratio_nonsense_syn_nHet = nHet_nonsense / nHet_synonymous,
          ratio_nonsense_syn_nNonRefHom = nNonRefHom_nonsense / nNonRefHom_synonymous,
          `Replicate line` = line_rep,
          Treatment = treatment) %>% 
  select(sample, treatment, `Replicate line`, sex, ratio_nonsense_syn_totalsnps, ratio_nonsense_syn_nHet, ratio_nonsense_syn_nNonRefHom) %>% 
  pivot_longer(cols=c(ratio_nonsense_syn_totalsnps, ratio_nonsense_syn_nNonRefHom, ratio_nonsense_syn_nHet), names_to = "zygosity") %>% 
  rename(`SS regime` = treatment)

dd_snpeff_X <- bind_rows(dd_nonsense_autosomes, dd_nonsense_x, dd_synonymous_autosomes,
                         dd_synonymous_x) %>% 
  left_join(population_info) %>% 
  filter(snpEff_impact %in% c("nonsense", "synonymous")) %>% 
  filter(chrom=="X") %>% 
  mutate(total_snps = nHet + nNonRefHom) %>% 
  mutate(snpEff_impact = case_when( snpEff_impact=="nonsense" ~ "nonsense",
                                    TRUE ~ snpEff_impact),
         treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Polyandry" ))) %>% 
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, snpEff_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=snpEff_impact, values_from=total_snps) %>% 
  mutate( ratio_nonsense_syn_totalsnps = nonsense / synonymous,
          `Replicate line` = line_rep,
          Treatment = treatment)


bind_rows(dd_nonsense_autosomes, dd_nonsense_x, dd_synonymous_autosomes,
          dd_synonymous_x) %>% 
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
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(-chrom, -nRefHom, -seq_run, -treat_rep) %>% 
  pivot_wider(names_from=snpEff_impact, values_from=c(total_snps, nHet, nNonRefHom)) %>% 
  mutate( ratio_nonsense_syn_totalsnps = total_snps_nonsense / total_snps_synonymous,
          ratio_nonsense_syn_nHet = nHet_nonsense / nHet_synonymous,
          ratio_nonsense_syn_nNonRefHom = nNonRefHom_nonsense / nNonRefHom_synonymous,
          `Replicate line` = line_rep,
          Treatment = treatment) %>% 
  select(sample, treatment, line_rep, sex, ratio_nonsense_syn_totalsnps, ratio_nonsense_syn_nHet, ratio_nonsense_syn_nNonRefHom) %>% 
  group_by(treatment, line_rep) %>% 
  summarise(mean_ratio_non_syn_totalsnps = mean(ratio_nonsense_syn_totalsnps),
            mean_ratio_non_syn_nHet = mean(ratio_nonsense_syn_nHet),
            mean_ratio_non_syn_nNonRefHom = mean(ratio_nonsense_syn_nNonRefHom)) %>% 
  write_rds(file="/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/lumley2015_godwin2020/inbreeding_resiliance_old_data/non_syn_ratios_by_line.rds",
  )



plot_nonsense_autosomes <- dd_snpeff_autosomes %>% 
  filter(zygosity == "ratio_nonsense_syn_totalsnps") %>% 
  ggplot(aes(x=.data[["SS regime"]], y=value, fill=.data[["SS regime"]])) +
  geom_jitter(aes(shape=.data[["Replicate line"]]),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Nonsense",
       y="",
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0.0034, 0.00490))

plot_nonsense_autosomes_BYZYG <- dd_snpeff_autosomes %>% 
  filter(zygosity != "ratio_nonsense_syn_totalsnps") %>% 
  mutate(zygosity = case_when(zygosity == "ratio_nonsense_syn_nHet" ~ "Het",
                              zygosity == "ratio_nonsense_syn_nNonRefHom" ~ "Hom")) %>% 
  ggplot(aes(x=`SS regime`, y=value, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Nonsense",
       y=expression(N[X]~" / Syn ratio"),
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~zygosity)


plot_nonsense_X <- bind_rows(dd_nonsense_autosomes, dd_nonsense_x, dd_synonymous_autosomes,
                             dd_synonymous_x) %>% 
  left_join(population_info) %>% 
  filter(snpEff_impact %in% c("nonsense", "synonymous")) %>% 
  filter(chrom=="X") %>% 
  mutate(total_snps = nHet + nNonRefHom) %>% 
  mutate(snpEff_impact = case_when( snpEff_impact=="nonsense" ~ "nonsense",
                                    TRUE ~ snpEff_impact),
         treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Polyandry" ))) %>% 
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, snpEff_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=snpEff_impact, values_from=total_snps) %>% 
  mutate( ratio_nonsense_syn_totalsnps = nonsense / synonymous,
          `Replicate line` = line_rep,
          `SS regime` = treatment) %>% 
  ggplot(aes(x=`SS regime`, y=ratio_nonsense_syn_totalsnps, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Nonsense",
       y="",
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0.001, 0.005))





plot_nonsense_autosomes_MvF <- bind_rows(dd_nonsense_autosomes, dd_nonsense_x, dd_synonymous_autosomes,
                                         dd_synonymous_x) %>% 
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
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, snpEff_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=snpEff_impact, values_from=total_snps) %>% 
  mutate( ratio_nonsense_syn_totalsnps = nonsense / synonymous,
          `Replicate line` = line_rep,
          `SS regime` = treatment,
          sex = case_when( sex == "male" ~ "Male",
                           sex == "female" ~ "Female")) %>% 
  ggplot(aes(x=`SS regime`, y=ratio_nonsense_syn_totalsnps, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Nonsense",
       y="",
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~sex)




plot_nonsense_X_MvF <- bind_rows(dd_nonsense_autosomes, dd_nonsense_x, dd_synonymous_autosomes,
                                 dd_synonymous_x) %>% 
  left_join(population_info) %>% 
  filter(snpEff_impact %in% c("nonsense", "synonymous")) %>% 
  filter(chrom=="X") %>% 
  mutate(total_snps = nHet + nNonRefHom) %>% 
  mutate(snpEff_impact = case_when( snpEff_impact=="nonsense" ~ "nonsense",
                                    TRUE ~ snpEff_impact),
         treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Polyandry" ))) %>% 
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, snpEff_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=snpEff_impact, values_from=total_snps) %>% 
  mutate( ratio_nonsense_syn_totalsnps = nonsense / synonymous,
          `Replicate line` = line_rep,
          `SS regime` = treatment,
          sex = case_when( sex == "male" ~ "Male",
                           sex == "female" ~ "Female")) %>% 
  ggplot(aes(x=`SS regime`, y=ratio_nonsense_syn_totalsnps, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Nonsense",
       y="",
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~sex)












#######################     GET THE RELEVANT SIFT DATA & PLOT - MISSENSE DELETERIOUS / TOLERATED
#######################


setwd("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/SIFT")

##################################  Get the info on all the samples
source("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/sample_info/generate_population.list.full.R")
population_info <- generate_ss_population.info.full.list("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/sample_info/sexual_selection_samplenames") %>%
  rename(sample = name) %>%
  mutate( line_rep = case_when( sample== "polym4" ~ "C",
                                TRUE ~ line_rep)) %>%
  mutate(treat_rep = paste(treatment, line_rep, sep="_"))





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



########### PULL DATA TOGETHER
##################

dd_sift_autosomes <- bind_rows(
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
  dd_missense_tol_autosomes,
  dd_missense_tol_X,
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
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_snps, nHet, nNonRefHom) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=c(total_snps, nHet, nNonRefHom)) %>% 
  mutate( ratio_misdel_syn_totalsnps = total_snps_deleterious / total_snps_synonymous,
          ratio_mistol_syn_totalsnps = total_snps_tolerated / total_snps_synonymous,
          ratio_misdel_syn_nHet = nHet_deleterious / nHet_synonymous,
          ratio_mistol_syn_nHet = nHet_tolerated / nHet_synonymous,
          ratio_misdel_syn_nNonRefHom = nNonRefHom_deleterious / nNonRefHom_synonymous,
          ratio_mistol_syn_nNonRefHom = nNonRefHom_tolerated / nNonRefHom_synonymous,
          `Replicate line` = line_rep,
          Treatment = treatment)

dd_sift_X <- bind_rows(
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
  dd_missense_tol_autosomes,
  dd_missense_tol_X,
  dd_synonymous_autosomes,
  dd_synonymous_X
) %>%
  left_join(population_info)  %>%
  filter(dataset %in% c("missense","synonymous"),
         SIFT_impact %in% c("deleterious", "tolerated", "synonymous")) %>%
  filter(chrom=="X") %>%
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


bind_rows(
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
  dd_missense_tol_autosomes,
  dd_missense_tol_X,
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
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_snps, nHet, nNonRefHom) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=c(total_snps, nHet, nNonRefHom)) %>% 
  mutate( ratio_misdel_syn_totalsnps = total_snps_deleterious / total_snps_synonymous,
          ratio_mistol_syn_totalsnps = total_snps_tolerated / total_snps_synonymous,
          ratio_misdel_syn_nHet = nHet_deleterious / nHet_synonymous,
          ratio_mistol_syn_nHet = nHet_tolerated / nHet_synonymous,
          ratio_misdel_syn_nNonRefHom = nNonRefHom_deleterious / nNonRefHom_synonymous,
          ratio_mistol_syn_nNonRefHom = nNonRefHom_tolerated / nNonRefHom_synonymous,
          Treatment = treatment) %>% 
  select(sample, treatment, line_rep, sex, ratio_misdel_syn_totalsnps, ratio_mistol_syn_totalsnps, ratio_misdel_syn_nHet, ratio_mistol_syn_nHet, ratio_misdel_syn_nNonRefHom, ratio_mistol_syn_nNonRefHom) %>% 
  group_by(treatment, line_rep) %>% 
  summarise( mean_ratio_misdel_syn_totalsnps = mean(ratio_misdel_syn_totalsnps),
             mean_ratio_mistol_syn_totalsnps = mean(ratio_mistol_syn_totalsnps),
             mean_ratio_misdel_syn_nHet = mean(ratio_misdel_syn_nHet),
             mean_ratio_mistol_syn_nHet = mean(ratio_mistol_syn_nHet),
             mean_ratio_misdel_syn_nNonRefHom = mean(ratio_misdel_syn_nNonRefHom),
             mean_ratio_mistol_syn_nNonRefHom = mean(ratio_mistol_syn_nNonRefHom)) %>% 
  write_rds(file="/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/lumley2015_godwin2020/inbreeding_resiliance_old_data/misdel_syn_ratios_by_line.rds",
  )



plot_misdel_autosomes <- bind_rows(
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
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
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=total_snps) %>% 
  mutate( ratio_misdel_syn_totalsnps = deleterious / synonymous,
          ratio_mistol_syn_totalsnps = tolerated / synonymous,
          `Replicate line` = line_rep,
          `SS regime` = treatment) %>% 
  ggplot(aes(x=`SS regime`, y=ratio_misdel_syn_totalsnps, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Missense deleterious ",
       y="",
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0.027, 0.0312))

plot_mistol_autosomes <- bind_rows(
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
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
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=total_snps) %>% 
  mutate( ratio_misdel_syn_totalsnps = deleterious / synonymous,
          ratio_mistol_syn_totalsnps = tolerated / synonymous,
          `Replicate line` = line_rep,
          `SS regime` = treatment) %>% 
  ggplot(aes(x=`SS regime`, y=ratio_mistol_syn_totalsnps, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Missense tolerated ",
       y=expression(N[X]~" / Syn ratio"),
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0.233, 0.248))


plot_misdel_autosomes_MvF <- bind_rows(
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
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
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=total_snps) %>% 
  mutate( ratio_misdel_syn_totalsnps = deleterious / synonymous,
          ratio_mistol_syn_totalsnps = tolerated / synonymous,
          `Replicate line` = line_rep,
          `SS regime` = treatment,
          sex = case_when( sex =="male" ~ "Male",
                           sex == "female" ~ "Female")) %>% 
  ggplot(aes(x=`SS regime`, y=ratio_misdel_syn_totalsnps, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Missense deleterious ",
       y="",
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~sex)


plot_mistol_autosomes_MvF <- bind_rows(
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
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
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=total_snps) %>% 
  mutate( ratio_misdel_syn_totalsnps = deleterious / synonymous,
          ratio_mistol_syn_totalsnps = tolerated / synonymous,
          `Replicate line` = line_rep,
          `SS regime` = treatment,
          sex = case_when( sex == "male" ~ "Male",
                           sex == "female" ~ "Female")) %>% 
  ggplot(aes(x=`SS regime`, y=ratio_mistol_syn_totalsnps, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Missense tolerated ",
       y=expression(N[X]~" / Syn ratio"),
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~sex)




############## PLOTTING X

plot_misdel_X <- bind_rows(
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
  dd_missense_tol_autosomes,
  dd_synonymous_autosomes,
  dd_synonymous_X
) %>% 
  left_join(population_info)  %>%  
  filter(dataset %in% c("missense","synonymous"),
         SIFT_impact %in% c("deleterious", "tolerated", "synonymous")) %>% 
  filter(chrom=="X") %>%
  filter(conf %in% c("high", "tol","syn")) %>% 
  mutate(total_snps = nHet + nNonRefHom) %>%
  mutate(treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Polyandry" ))) %>% 
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=total_snps) %>% 
  mutate( ratio_misdel_syn_totalsnps = deleterious / synonymous,
          `Replicate line` = line_rep,
          `SS regime` = treatment) %>% 
  ggplot(aes(x=`SS regime`, y=ratio_misdel_syn_totalsnps, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Missense deleterious ",
       y="",
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0.018, 0.0295))




plot_mistol_X <- bind_rows(
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
  dd_missense_tol_autosomes,
  dd_missense_tol_X,
  dd_synonymous_autosomes,
  dd_synonymous_X
) %>% 
  left_join(population_info)  %>%  
  filter(dataset %in% c("missense","synonymous"),
         SIFT_impact %in% c("deleterious", "tolerated", "synonymous")) %>% 
  filter(chrom=="X") %>%
  filter(conf %in% c("high", "tol","syn")) %>% 
  mutate(total_snps = nHet + nNonRefHom) %>%
  mutate(treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Polyandry" ))) %>% 
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=total_snps) %>% 
  mutate( ratio_mistol_syn_totalsnps = tolerated / synonymous,
          `Replicate line` = line_rep,
          `SS regime` = treatment) %>% 
  ggplot(aes(x=`SS regime`, y=ratio_mistol_syn_totalsnps, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5, outlier.shape = NA)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Missense tolerated ",
       y=expression(N[X]~" / Syn ratio"),
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5)) +
  coord_cartesian(ylim = c(0.18, 0.245))



plot_misdel_X_MvF <- bind_rows(
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
  dd_missense_tol_autosomes,
  dd_synonymous_autosomes,
  dd_synonymous_X
) %>% 
  left_join(population_info)  %>%  
  filter(dataset %in% c("missense","synonymous"),
         SIFT_impact %in% c("deleterious", "tolerated", "synonymous")) %>% 
  filter(chrom=="X") %>%
  filter(conf %in% c("high", "tol","syn")) %>% 
  mutate(total_snps = nHet + nNonRefHom) %>%
  mutate(treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Polyandry" ))) %>% 
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=total_snps) %>% 
  mutate( ratio_misdel_syn_totalsnps = deleterious / synonymous,
          `Replicate line` = line_rep,
          `SS regime` = treatment,
          sex = case_when( sex == "male" ~ "Male",
                           sex == "female" ~ "Female")) %>% 
  ggplot(aes(x=`SS regime`, y=ratio_misdel_syn_totalsnps, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Missense deleterious ",
       y="",
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~sex)


plot_mistol_X_MvF <- bind_rows(
  dd_missense_del_highconf_autosomes,
  dd_missense_del_highconf_X,
  dd_missense_tol_autosomes,
  dd_missense_tol_X,
  dd_synonymous_autosomes,
  dd_synonymous_X
) %>% 
  left_join(population_info)  %>%  
  filter(dataset %in% c("missense","synonymous"),
         SIFT_impact %in% c("deleterious", "tolerated", "synonymous")) %>% 
  filter(chrom=="X") %>%
  filter(conf %in% c("high", "tol","syn")) %>% 
  mutate(total_snps = nHet + nNonRefHom) %>%
  mutate(treatment = case_when( treatment=="GA1" ~ "Ancestral",
                                treatment=="f_biased" ~ "Female biased",
                                treatment=="m_biased" ~ "Male biased",
                                treatment=="mono" ~ "Monandry",
                                treatment=="poly" ~ "Polyandry"),
         treatment = factor(treatment, levels=c("Ancestral","Female biased","Male biased","Monandry","Polyandry" ))) %>% 
  filter(treatment %in% c("Monandry", "Polyandry")) %>% 
  select(sample, SIFT_impact, treatment, line_rep, sex, total_snps) %>% 
  pivot_wider(names_from=SIFT_impact, values_from=total_snps) %>% 
  mutate( ratio_mistol_syn_totalsnps = tolerated / synonymous,
          `Replicate line` = line_rep,
          `SS regime` = treatment,
          sex = case_when( sex == "male" ~ "Male",
                           sex == "female" ~ "Female")) %>% 
  ggplot(aes(x=`SS regime`, y=ratio_mistol_syn_totalsnps, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), alpha = 0.4) +
  geom_boxplot(alpha=0.5)+
  theme_bw() +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(title="Missense tolerated ",
       y=expression(N[X]~" / Syn ratio"),
       x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~sex)



plots_autosomes <- ggarrange(plot_mistol_autosomes, plot_misdel_autosomes, plot_nonsense_autosomes, nrow=1, common.legend=T, legend="bottom")
ggsave("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/varCat_syn_ratio_SNPS_autosomes.png",
       plots_autosomes,
       device="png",
       units="mm",
       width=200,
       height=100)

ggsave("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/varCat_syn_ratio_SNPS_autosomes_byZygosity.png",
       plot_nonsense_autosomes_BYZYG,
       device="png",
       units="mm",
       width=120,
       height=100)

plots_X <- ggarrange(plot_mistol_X, plot_misdel_X, plot_nonsense_X, nrow=1, common.legend=T, legend="bottom")
ggsave("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/varCat_syn_ratio_SNPS_X.png",
       plots_X,
       device="png",
       units="mm",
       width=200,
       height=100)


plot_autosomes_by_sex <- ggarrange(plot_mistol_autosomes_MvF,
                                   plot_misdel_autosomes_MvF,
                                   plot_nonsense_autosomes_MvF,
                                   nrow=1, common.legend = T, legend="bottom")
ggsave("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/varCat_syn_ratio_SNPS_autosomes_by_sex.png",
       plot_autosomes_by_sex,
       device="png",
       units="mm",
       width=200,
       height=100)

plot_X_by_sex <- ggarrange(plot_mistol_X_MvF,
                           plot_misdel_X_MvF,
                           plot_nonsense_X_MvF,
                           nrow=1, common.legend = T, legend="bottom")
ggsave("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/varCat_syn_ratio_SNPS_X_by_sex.png",
       plot_X_by_sex,
       device="png",
       units="mm",
       width=200,
       height=100)



##################################
##################################
################################## MODELLING

# ###### check distributions
# 
# dd_snpeff_autosomes %>%
#   ggplot(aes(ratio_nonsense_syn_totalsnps)) +
#   geom_histogram() +
#   theme_bw() +
#   facet_wrap(~treatment)
# 
# dd_sift_autosomes %>%
#   ggplot(aes(ratio_misdel_syn_totalsnps)) +
#   geom_histogram() +
#   theme_bw() +
# facet_wrap(~treatment)
# 
# dd_sift_autosomes %>%
#   ggplot(aes(ratio_mistol_syn_totalsnps)) +
#   geom_histogram() +
#   theme_bw() +
# facet_wrap(~treatment)
# 
# 
# dd_snpeff_X %>%
#   ggplot(aes(ratio_nonsense_syn_totalsnps)) +
#   geom_histogram() +
#   theme_bw() +
# facet_wrap(~treatment)
# 
# dd_sift_X %>%
#   ggplot(aes(ratio_misdel_syn_totalsnps)) +
#   geom_histogram() +
#   theme_bw() +
# facet_wrap(~treatment)
# 
# dd_sift_X %>%
#   ggplot(aes(ratio_mistol_syn_totalsnps)) +
#   geom_histogram() +
#   theme_bw() +
# facet_wrap(~treatment)


##################. AUTOSOMES

################# nonsense
dd_snpeff_autosomes <- dd_snpeff_autosomes %>%
  filter(treatment %in% c("Monandry", "Polyandry"))
m_non_auto_int <- lmer(ratio_nonsense_syn_totalsnps ~ treatment + sex + treatment*sex +(1|line_rep), data=dd_snpeff_autosomes)
summary(m_non_auto_int)

m_non_auto <- lmer(ratio_nonsense_syn_totalsnps ~ treatment + sex +(1|line_rep), data=dd_snpeff_autosomes)
summary(m_non_auto)

# emm_non <- emmeans(m_non_auto, pairwise ~ treatment)
# emm_non

################# missense deleterious
dd_sift_autosomes <- dd_sift_autosomes %>% 
  filter(treatment %in% c("Monandry", "Polyandry"))

m_misdel_auto_int <- lmer(ratio_misdel_syn_totalsnps ~ treatment + sex + treatment*sex +(1|line_rep), data=dd_sift_autosomes)
summary(m_misdel_auto_int)

m_misdel_auto <- lmer(ratio_misdel_syn_totalsnps ~ treatment + sex +(1|line_rep), data=dd_sift_autosomes)
summary(m_misdel_auto)

# emm_misdel <- emmeans(m_misdel_auto, pairwise ~ treatment)
# emm_misdel


################# missense tolerated

m_mistol_auto_int <- lmer(ratio_mistol_syn_totalsnps ~ treatment + sex + treatment*sex +(1|line_rep), data=dd_sift_autosomes)
summary(m_mistol_auto_int)

m_mistol_auto <- lmer(ratio_mistol_syn_totalsnps ~ treatment + sex +(1|line_rep), data=dd_sift_autosomes)
summary(m_mistol_auto)

# emm_mistol <- emmeans(m_mistol_auto, pairwise ~ treatment)
# emm_mistol






#################  X CHROM

################# nonsense
dd_snpeff_X <- dd_snpeff_X %>% 
  filter( treatment %in% c("Monandry","Polyandry"))
m_non_X_int <- lmer(ratio_nonsense_syn_totalsnps ~ treatment + sex + treatment*sex +(1|line_rep), data=dd_snpeff_X)
summary(m_non_X_int)

m_non_X <- lmer(ratio_nonsense_syn_totalsnps ~ treatment + sex +(1|line_rep), data=dd_snpeff_X)
summary(m_non_X)

# emm_non_X <- emmeans(m_non_X, pairwise ~ treatment)
# emm_non_X



################# missense deleterious
dd_sift_X <- dd_sift_X %>% 
  filter( treatment %in% c("Monandry","Polyandry"))

m_misdel_X_int <- lmer(ratio_misdel_syn_totalsnps ~ treatment + sex + treatment*sex +(1|line_rep), data=dd_sift_X)
summary(m_misdel_X_int)

m_misdel_X <- lmer(ratio_misdel_syn_totalsnps ~ treatment + sex +(1|line_rep), data=dd_sift_X)
summary(m_misdel_X)

# emm_misdel_X <- emmeans(m_misdel_X, pairwise ~ treatment)
# emm_misdel_X


################# missense tolerated
m_mistol_X_int <- lmer(ratio_mistol_syn_totalsnps ~ treatment + sex + treatment*sex +(1|line_rep), data=dd_sift_X)
summary(m_mistol_X_int)

m_mistol_X <- lmer(ratio_mistol_syn_totalsnps ~ treatment + sex +(1|line_rep), data=dd_sift_X)
summary(m_mistol_X)

# emm_mistol_X <- emmeans(m_mistol_X, pairwise ~ treatment)
# emm_mistol_X
