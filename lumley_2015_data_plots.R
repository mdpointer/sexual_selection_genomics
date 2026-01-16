library(tidyverse)
library(survminer)
library(coxme)


## data from new SIFT analysis, ratio nonsense / synonymous variants
dd_snpeff <- read_rds("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/lumley2015_godwin2020/inbreeding_resiliance_old_data/non_syn_ratios_by_line.rds") %>% 
  mutate(key_field = paste(treatment, line_rep, sep="_"))

dd_SIFT <- read_rds("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/lumley2015_godwin2020/inbreeding_resiliance_old_data/misdel_syn_ratios_by_line.rds") %>% 
  mutate(key_field = paste(treatment, line_rep, sep="_"))


# selection B experiment (mono / poly)
# there are 8 reps per line
dd <- read.csv("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/lumley2015_godwin2020/inbreeding_resiliance_old_data/lumley_selB_survival_data.csv") %>% 
  pivot_longer( F1:F20, names_to = "gen", values_to = "surviving") %>% 
  mutate(gen = as.double(str_sub(gen, 2, length(gen))),
         rep = paste(Treatment, Line, sep="_")) %>% 
  mutate(key_field = paste(Treatment, Line, sep="_"))


dd <- full_join(dd, dd_snpeff, by="key_field") %>% 
  left_join(dd_SIFT, by="key_field") %>% 
  select(-key_field, -treatment.x, -line_rep.x, -treatment.y, -line_rep.y)


pop_extinctions <- dd %>%
  mutate(pop_alive = round(surviving * 8)) %>%
  arrange(rep, gen) %>%
  group_by(rep) %>%
  # Step 2: Count how many populations went extinct per generation
  mutate(extinct_this_gen = dplyr::lag(pop_alive, default = 8) - pop_alive) %>%
  ungroup() %>%
  # Step 3: Filter to generations where at least one extinction occurred
  filter(extinct_this_gen > 0) %>%
  # Step 4: Repeat the generation number for each population that died
  uncount(extinct_this_gen, .id = "pop_id") %>%
  dplyr::rename(extinction_gen = gen) %>%
  mutate(status = 1)

# Step 5: Add censored populations (those that survived to end)
max_gen <- max(dd$gen)
## No popultions are censored, so this produces an empty dataframe that we don't use further
censored <- pop_extinctions %>%
  count(rep, name = "extinct_n") %>%
  mutate(survived_n = 8 - extinct_n) %>%
  filter(survived_n > 0) %>%
  uncount(survived_n) %>%
  mutate(
    extinction_gen = max_gen,
    status = 0
  )

# Step 6: Combine extinct and censored populations
survival_data <- pop_extinctions %>%
  left_join(dd %>% select(rep, Treatment, mean_ratio_non_syn_totalsnps, mean_ratio_non_syn_nHet, mean_ratio_non_syn_nNonRefHom, mean_ratio_mistol_syn_totalsnps, mean_ratio_mistol_syn_nHet, mean_ratio_mistol_syn_nNonRefHom, mean_ratio_misdel_syn_totalsnps, mean_ratio_misdel_syn_nHet, mean_ratio_misdel_syn_nNonRefHom) %>% distinct(), by = "rep") %>% 
  select(-mean_ratio_non_syn_totalsnps.y, -mean_ratio_non_syn_nHet.y, -mean_ratio_non_syn_nNonRefHom.y, -mean_ratio_mistol_syn_totalsnps.y, -mean_ratio_mistol_syn_nHet.y, -mean_ratio_mistol_syn_nNonRefHom.y, -mean_ratio_misdel_syn_totalsnps.y, -mean_ratio_misdel_syn_nHet.y, -mean_ratio_misdel_syn_nNonRefHom.y,-Treatment.y) %>% 
  dplyr::rename(mean_ratio_non_syn_totalsnps = mean_ratio_non_syn_totalsnps.x,
         mean_ratio_non_syn_nHet = mean_ratio_non_syn_nHet.x,
         mean_ratio_non_syn_nNonRefHom = mean_ratio_non_syn_nNonRefHom.x,
         mean_ratio_mistol_syn_totalsnps = mean_ratio_mistol_syn_totalsnps.x, 
         mean_ratio_mistol_syn_nHet = mean_ratio_mistol_syn_nHet.x, 
         mean_ratio_mistol_syn_nNonRefHom = mean_ratio_mistol_syn_nNonRefHom.x, 
         mean_ratio_misdel_syn_totalsnps = mean_ratio_misdel_syn_totalsnps.x, 
         mean_ratio_misdel_syn_nHet = mean_ratio_misdel_syn_nHet.x, 
         mean_ratio_misdel_syn_nNonRefHom = mean_ratio_misdel_syn_nNonRefHom.x,
         Treatment = Treatment.x)

## mutate the cols to use in model, considering different contributions of nHet and nNonRefHom
survival_data <- survival_data %>% 
  mutate( 
    ## no weighting
          load_rec_noWeight = mean_ratio_misdel_syn_nNonRefHom + mean_ratio_non_syn_nNonRefHom,          # variants are recessive, no contribution from Hets
          load_parrec_noWeight = 0.25*(mean_ratio_misdel_syn_nHet + mean_ratio_non_syn_nHet) +             # Hets contribute 0.25
            mean_ratio_misdel_syn_nNonRefHom + 
            mean_ratio_non_syn_nNonRefHom,
          load_codom_noWeight = 0.5*(mean_ratio_misdel_syn_nHet + mean_ratio_non_syn_nHet) +             # Hets contribute 0.5
            mean_ratio_misdel_syn_nNonRefHom + 
            mean_ratio_non_syn_nNonRefHom,
   ## nonsense 2* deleteriousness of missense
          load_rec_2Weight = mean_ratio_misdel_syn_nNonRefHom + 2*mean_ratio_non_syn_nNonRefHom,       # recessive, plus nonsense twice as deleterious as missense
   load_parrec_2Weight = 0.25*(mean_ratio_misdel_syn_nHet + 2*mean_ratio_non_syn_nHet) +            # Hets contribute half & nonsense twice as del as missense
     mean_ratio_misdel_syn_nNonRefHom + 
     2*mean_ratio_non_syn_nNonRefHom,
          load_codom_2Weight = 0.5*(mean_ratio_misdel_syn_nHet + 2*mean_ratio_non_syn_nHet) +            # Hets contribute half & nonsense twice as del as missense
            mean_ratio_misdel_syn_nNonRefHom + 
            2*mean_ratio_non_syn_nNonRefHom,
   ## nonsense 5* deleteriousness of missense
          load_rec_5Weight = mean_ratio_misdel_syn_nNonRefHom + 5*mean_ratio_non_syn_nNonRefHom,       # recessive, plus nonsense 5* as deleterious as missense
   load_parrec_5Weight = 0.5*(mean_ratio_misdel_syn_nHet + 5*mean_ratio_non_syn_nHet) +            # Hets contribute half & nonsense 5* as del as missense
     mean_ratio_misdel_syn_nNonRefHom + 
     5*mean_ratio_non_syn_nNonRefHom,
          load_codom_5Weight = 0.5*(mean_ratio_misdel_syn_nHet + 5*mean_ratio_non_syn_nHet) +            # Hets contribute half & nonsense 5* as del as missense
            mean_ratio_misdel_syn_nNonRefHom + 
            5*mean_ratio_non_syn_nNonRefHom
  )





########
#######      MODELLING IN A PROP HAZARDS FRAMEWORK - CAN ADD REPLICATE AS A RANDOM FACTOR

###### MODELLING THE EFFECT OF TREATMENT VS THE EFFECT OF LOAD, UNDER A FEW DIFFERENT METHODS OF QUANTIFYING LOAD:


# ###### VERIFYING THAT THE ASSUMPTION OF PROPORTIONAL HAZARDS IS SATISFIED FOR THESE DATA - LOOKS GOOD
# cox_basic <- coxph(Surv(extinction_gen, status) ~ Treatment + load_rec_noWeight, data = survival_data)
# ph_test <- cox.zph(cox_basic)
# print(ph_test)
# 
# coxme_model <- coxme(Surv(extinction_gen, status) ~ Treatment + load_rec_noWeight + (1 | rep), data = survival_data)
# summary(coxme_model)
# 
# cox_cluster <- coxph(Surv(extinction_gen, status) ~ Treatment + load_rec_noWeight + cluster(rep), data = survival_data)
# cox.zph(cox_cluster)



###### 1) TOTAL DOMINANCE OF REF ALLELE (ONLY NON REF HOMS CONSIDERED) & NONSENSE AND MISSENSE GIVEN EQUAL WEIGHT
# Model 1: Treatment only
model_treatment <- coxme(Surv(extinction_gen, status) ~ Treatment + (1|rep),
                           data = survival_data)

# #######Test that the proportional hazards assumption is not violated
# ph_model <- coxph(Surv(extinction_gen, status) ~ Treatment, data = survival_data)
# # Test proportional hazards assumption
# ph_test <- cox.zph(ph_model)
# ph_test


summary(model_treatment)
# Model 2: Mutation load nonsense only
model_load_rec_noWeight <- coxme(Surv(extinction_gen, status) ~ load_rec_noWeight + (1|rep),
                                   data = survival_data)
summary(model_load_rec_noWeight)

# Model 3: Treatment + Mutation load
model_load_rec_noWeight_both <- coxme(Surv(extinction_gen, status) ~ Treatment + load_rec_noWeight + (1|rep),
                                        data = survival_data)
summary(model_load_rec_noWeight_both)

# Extract log-likelihoods using logLik()
loglik_model_treatment <- as.numeric(logLik(model_treatment))
loglik_model_load_rec_noWeight <- as.numeric(logLik(model_load_rec_noWeight))
loglik_model_load_rec_noWeight_both <- as.numeric(logLik(model_load_rec_noWeight_both))

# Count number of fixed effects parameters
n_fixed_treatment <- length(fixef(model_treatment))
n_fixed_load <- length(fixef(model_load_rec_noWeight))
n_fixed_both <- length(fixef(model_load_rec_noWeight_both))

# Count number of random effects variance parameters
# coxme stores variance components in model$variance (named vector)
n_random_treatment <- length(model_treatment$variance)
n_random_load <- length(model_load_rec_noWeight$variance)
n_random_both <- length(model_load_rec_noWeight_both$variance)

# Total number of parameters = fixed effects + random effects variances
npar_model_treatment <- n_fixed_treatment + n_random_treatment
npar_model_load_rec_noWeight <- n_fixed_load + n_random_load
npar_model_load_rec_noWeight_both <- n_fixed_both + n_random_both

# Calculate AIC manually
aic_model_treatment <- -2 * loglik_model_treatment + 2 * npar_model_treatment
aic_model_load_rec_noWeight <- -2 * loglik_model_load_rec_noWeight + 2 * npar_model_load_rec_noWeight
aic_model_load_rec_noWeight_both <- -2 * loglik_model_load_rec_noWeight_both + 2 * npar_model_load_rec_noWeight_both

# Make a comparison table
model_comparison_load_rec_noWeight <- data.frame(
  Load_scenario = "load_rec_noWeight",
  Model = c("Treatment", "Mutation Load", "Treatment + Load"),
  LogLik = c(loglik_model_treatment, loglik_model_load_rec_noWeight, loglik_model_load_rec_noWeight_both),
  n_params = c(npar_model_treatment, npar_model_load_rec_noWeight, npar_model_load_rec_noWeight_both),
  AIC = c(aic_model_treatment, aic_model_load_rec_noWeight, aic_model_load_rec_noWeight_both)
)

print(model_comparison_load_rec_noWeight)





###### 2) PARTIAL RECESSIVENESS (NON REF HOMS COUNT DOUBLE HET) & MISSENSE AND NONSENSE GIVEN EQUAL WEIGHT
# Model 1: Treatment only
model_treatment <- coxme(Surv(extinction_gen, status) ~ Treatment + (1|rep),
                         data = survival_data)
summary(model_treatment)
# Model 2: Mutation load nonsense only
model_load_parrec_noWeight <- coxme(Surv(extinction_gen, status) ~ load_parrec_noWeight + (1|rep),
                                   data = survival_data)
summary(model_load_parrec_noWeight)
# Model 3: Treatment + Mutation load
model_load_parrec_noWeight_both <- coxme(Surv(extinction_gen, status) ~ Treatment + load_parrec_noWeight + (1|rep),
                                        data = survival_data)
summary(model_load_parrec_noWeight_both)
# Extract log-likelihoods using logLik()
loglik_model_treatment <- as.numeric(logLik(model_treatment))
loglik_model_load_parrec_noWeight <- as.numeric(logLik(model_load_parrec_noWeight))
loglik_model_load_parrec_noWeight_both <- as.numeric(logLik(model_load_parrec_noWeight_both))

# Count number of fixed effects parameters
n_fixed_treatment <- length(fixef(model_treatment))
n_fixed_load <- length(fixef(model_load_parrec_noWeight))
n_fixed_both <- length(fixef(model_load_parrec_noWeight_both))

# Count number of random effects variance parameters
n_random_treatment <- length(model_treatment$variance)
n_random_load <- length(model_load_parrec_noWeight$variance)
n_random_both <- length(model_load_parrec_noWeight_both$variance)

# Total number of parameters = fixed effects + random effects variances
npar_model_treatment <- n_fixed_treatment + n_random_treatment
npar_model_load_parrec_noWeight <- n_fixed_load + n_random_load
npar_model_load_parrec_noWeight_both <- n_fixed_both + n_random_both

# Calculate AIC manually
aic_model_treatment <- -2 * loglik_model_treatment + 2 * npar_model_treatment
aic_model_load_parrec_noWeight <- -2 * loglik_model_load_parrec_noWeight + 2 * npar_model_load_parrec_noWeight
aic_model_load_parrec_noWeight_both <- -2 * loglik_model_load_parrec_noWeight_both + 2 * npar_model_load_parrec_noWeight_both

# Make a comparison table
model_comparison_load_parrec_noWeight <- data.frame(
  Load_scenario = "load_parrec_noWeight",
  Model = c("Treatment", "Mutation Load", "Treatment + Load"),
  LogLik = c(loglik_model_treatment, loglik_model_load_parrec_noWeight, loglik_model_load_parrec_noWeight_both),
  n_params = c(npar_model_treatment, npar_model_load_parrec_noWeight, npar_model_load_parrec_noWeight_both),
  AIC = c(aic_model_treatment, aic_model_load_parrec_noWeight, aic_model_load_parrec_noWeight_both)
)

print(model_comparison_load_parrec_noWeight)



###### 3) CO-DOMINANCE (NON REF HOMS COUNT DOUBLE HET) & MISSENSE AND NONSENSE GIVEN EQUAL WEIGHT
# Model 1: Treatment only
model_treatment <- coxme(Surv(extinction_gen, status) ~ Treatment + (1|rep),
                         data = survival_data)

# Model 2: Mutation load nonsense only
model_load_codom_noWeight <- coxme(Surv(extinction_gen, status) ~ load_codom_noWeight + (1|rep),
                                   data = survival_data)
summary(model_load_codom_noWeight)
# Model 3: Treatment + Mutation load
model_load_codom_noWeight_both <- coxme(Surv(extinction_gen, status) ~ Treatment + load_codom_noWeight + (1|rep),
                                        data = survival_data)
summary(model_load_codom_noWeight_both)
# Extract log-likelihoods using logLik()
loglik_model_treatment <- as.numeric(logLik(model_treatment))
loglik_model_load_codom_noWeight <- as.numeric(logLik(model_load_codom_noWeight))
loglik_model_load_codom_noWeight_both <- as.numeric(logLik(model_load_codom_noWeight_both))

# Count number of fixed effects parameters
n_fixed_treatment <- length(fixef(model_treatment))
n_fixed_load <- length(fixef(model_load_codom_noWeight))
n_fixed_both <- length(fixef(model_load_codom_noWeight_both))

# Count number of random effects variance parameters
n_random_treatment <- length(model_treatment$variance)
n_random_load <- length(model_load_codom_noWeight$variance)
n_random_both <- length(model_load_codom_noWeight_both$variance)

# Total number of parameters = fixed effects + random effects variances
npar_model_treatment <- n_fixed_treatment + n_random_treatment
npar_model_load_codom_noWeight <- n_fixed_load + n_random_load
npar_model_load_codom_noWeight_both <- n_fixed_both + n_random_both

# Calculate AIC manually
aic_model_treatment <- -2 * loglik_model_treatment + 2 * npar_model_treatment
aic_model_load_codom_noWeight <- -2 * loglik_model_load_codom_noWeight + 2 * npar_model_load_codom_noWeight
aic_model_load_codom_noWeight_both <- -2 * loglik_model_load_codom_noWeight_both + 2 * npar_model_load_codom_noWeight_both

# Make a comparison table
model_comparison_load_codom_noWeight <- data.frame(
  Load_scenario = "load_codom_noWeight",
  Model = c("Treatment", "Mutation Load", "Treatment + Load"),
  LogLik = c(loglik_model_treatment, loglik_model_load_codom_noWeight, loglik_model_load_codom_noWeight_both),
  n_params = c(npar_model_treatment, npar_model_load_codom_noWeight, npar_model_load_codom_noWeight_both),
  AIC = c(aic_model_treatment, aic_model_load_codom_noWeight, aic_model_load_codom_noWeight_both)
)

print(model_comparison_load_codom_noWeight)



###### 4) TOTAL DOMINANCE OF REF ALLELE (ONLY NON REF HOMS CONSIDERED) & NONSENSE COUNT 2X MISSENSE
# Model 1: Treatment only
model_treatment <- coxme(Surv(extinction_gen, status) ~ Treatment + (1|rep),
                         data = survival_data)

# Model 2: Mutation load nonsense only
model_load_rec_2Weight <- coxme(Surv(extinction_gen, status) ~ load_rec_2Weight + (1|rep),
                                data = survival_data)
summary(model_load_rec_2Weight)
# Model 3: Treatment + Mutation load
model_load_rec_2Weight_both <- coxme(Surv(extinction_gen, status) ~ Treatment + load_rec_2Weight + (1|rep),
                                     data = survival_data)
summary(model_load_rec_2Weight_both)
# Extract log-likelihoods using logLik()
loglik_model_treatment <- as.numeric(logLik(model_treatment))
loglik_model_load_rec_2Weight <- as.numeric(logLik(model_load_rec_2Weight))
loglik_model_load_rec_2Weight_both <- as.numeric(logLik(model_load_rec_2Weight_both))

# Count number of fixed effects parameters
n_fixed_treatment <- length(fixef(model_treatment))
n_fixed_load <- length(fixef(model_load_rec_2Weight))
n_fixed_both <- length(fixef(model_load_rec_2Weight_both))

# Count number of random effects variance parameters
n_random_treatment <- length(model_treatment$variance)
n_random_load <- length(model_load_rec_2Weight$variance)
n_random_both <- length(model_load_rec_2Weight_both$variance)

# Total number of parameters = fixed effects + random effects variances
npar_model_treatment <- n_fixed_treatment + n_random_treatment
npar_model_load_rec_2Weight <- n_fixed_load + n_random_load
npar_model_load_rec_2Weight_both <- n_fixed_both + n_random_both

# Calculate AIC manually
aic_model_treatment <- -2 * loglik_model_treatment + 2 * npar_model_treatment
aic_model_load_parrec_2Weight <- -2 * loglik_model_load_rec_2Weight + 2 * npar_model_load_rec_2Weight
aic_model_load_parrec_2Weight_both <- -2 * loglik_model_load_rec_2Weight_both + 2 * npar_model_load_rec_2Weight_both

# Make a comparison table
model_comparison_load_rec_2Weight <- data.frame(
  Load_scenario = "load_rec_2Weight",
  Model = c("Treatment", "Mutation Load", "Treatment + Load"),
  LogLik = c(loglik_model_treatment, loglik_model_load_rec_2Weight, loglik_model_load_rec_2Weight_both),
  n_params = c(npar_model_treatment, npar_model_load_rec_2Weight, npar_model_load_rec_2Weight_both),
  AIC = c(aic_model_treatment, aic_model_load_rec_2Weight, aic_model_load_rec_2Weight_both)
)

print(model_comparison_load_parrec_2Weight)



###### 5) PARTIAL RECESSIVENESS (NON REF HOMS COUNT DOUBLE HET)  & NONSENSE COUNT 2X MISSENSE
# Model 1: Treatment only
model_treatment <- coxme(Surv(extinction_gen, status) ~ Treatment + (1|rep),
                         data = survival_data)
summary(model_treatment)

# Model 2: Mutation load nonsense only
model_load_parrec_2Weight <- coxme(Surv(extinction_gen, status) ~ load_parrec_2Weight + (1|rep),
                                    data = survival_data)
summary(model_load_parrec_2Weight)

# Model 3: Treatment + Mutation load
model_load_parrec_2Weight_both <- coxme(Surv(extinction_gen, status) ~ Treatment + load_parrec_2Weight + (1|rep),
                                         data = survival_data)
summary(model_load_parrec_2Weight_both)

# Extract log-likelihoods using logLik()
loglik_model_treatment <- as.numeric(logLik(model_treatment))
loglik_model_load_parrec_2Weight <- as.numeric(logLik(model_load_parrec_2Weight))
loglik_model_load_parrec_2Weight_both <- as.numeric(logLik(model_load_parrec_2Weight_both))

# Count number of fixed effects parameters
n_fixed_treatment <- length(fixef(model_treatment))
n_fixed_load <- length(fixef(model_load_parrec_2Weight))
n_fixed_both <- length(fixef(model_load_parrec_2Weight_both))

# Count number of random effects variance parameters
n_random_treatment <- length(model_treatment$variance)
n_random_load <- length(model_load_parrec_2Weight$variance)
n_random_both <- length(model_load_parrec_2Weight_both$variance)

# Total number of parameters = fixed effects + random effects variances
npar_model_treatment <- n_fixed_treatment + n_random_treatment
npar_model_load_parrec_2Weight <- n_fixed_load + n_random_load
npar_model_load_parrec_2Weight_both <- n_fixed_both + n_random_both

# Calculate AIC manually
aic_model_treatment <- -2 * loglik_model_treatment + 2 * npar_model_treatment
aic_model_load_parrec_2Weight <- -2 * loglik_model_load_parrec_2Weight + 2 * npar_model_load_parrec_2Weight
aic_model_load_parrec_2Weight_both <- -2 * loglik_model_load_parrec_2Weight_both + 2 * npar_model_load_parrec_2Weight_both

# Make a comparison table
model_comparison_load_parrec_2Weight <- data.frame(
  Load_scenario = "load_parrec_2Weight",
  Model = c("Treatment", "Mutation Load", "Treatment + Load"),
  LogLik = c(loglik_model_treatment, loglik_model_load_parrec_2Weight, loglik_model_load_parrec_2Weight_both),
  n_params = c(npar_model_treatment, npar_model_load_parrec_2Weight, npar_model_load_parrec_2Weight_both),
  AIC = c(aic_model_treatment, aic_model_load_parrec_2Weight, aic_model_load_parrec_2Weight_both)
)

print(model_comparison_load_parrec_2Weight)





###### 6) CO-DOMINANCE (NON REF HOMS COUNT DOUBLE HET) & NONSENSE COUNT 2X MISSENSE
# Model 1: Treatment only
model_treatment <- coxme(Surv(extinction_gen, status) ~ Treatment + (1|rep),
                         data = survival_data)

# Model 2: Mutation load nonsense only
model_load_codom_2Weight <- coxme(Surv(extinction_gen, status) ~ load_codom_2Weight + (1|rep),
                                  data = survival_data)
summary(model_load_codom_2Weight)
# Model 3: Treatment + Mutation load
model_load_codom_2Weight_both <- coxme(Surv(extinction_gen, status) ~ Treatment + load_codom_2Weight + (1|rep),
                                       data = survival_data)
summary(model_load_codom_2Weight_both)
# Extract log-likelihoods
loglik_model_treatment <- as.numeric(logLik(model_treatment))
loglik_model_load_codom_2Weight <- as.numeric(logLik(model_load_codom_2Weight))
loglik_model_load_codom_2Weight_both <- as.numeric(logLik(model_load_codom_2Weight_both))

# Number of fixed effects parameters
n_fixed_treatment <- length(fixef(model_treatment))
n_fixed_load <- length(fixef(model_load_codom_2Weight))
n_fixed_both <- length(fixef(model_load_codom_2Weight_both))

# Number of random effects variance parameters
n_random_treatment <- length(model_treatment$variance)
n_random_load <- length(model_load_codom_2Weight$variance)
n_random_both <- length(model_load_codom_2Weight_both$variance)

# Total number of parameters = fixed + random effects variance parameters
npar_model_treatment <- n_fixed_treatment + n_random_treatment
npar_model_load_codom_2Weight <- n_fixed_load + n_random_load
npar_model_load_codom_2Weight_both <- n_fixed_both + n_random_both

# Calculate AIC manually
aic_model_treatment <- -2 * loglik_model_treatment + 2 * npar_model_treatment
aic_model_load_codom_2Weight <- -2 * loglik_model_load_codom_2Weight + 2 * npar_model_load_codom_2Weight
aic_model_load_codom_2Weight_both <- -2 * loglik_model_load_codom_2Weight_both + 2 * npar_model_load_codom_2Weight_both

# Create comparison table
model_comparison_load_codom_2Weight <- data.frame(
  Load_scenario = "load_codom_2Weight",
  Model = c("Treatment", "Mutation Load", "Treatment + Load"),
  LogLik = c(loglik_model_treatment, loglik_model_load_codom_2Weight, loglik_model_load_codom_2Weight_both),
  n_params = c(npar_model_treatment, npar_model_load_codom_2Weight, npar_model_load_codom_2Weight_both),
  AIC = c(aic_model_treatment, aic_model_load_codom_2Weight, aic_model_load_codom_2Weight_both)
)

print(model_comparison_load_codom_2Weight)



###### 7) TOTAL DOMINANCE OF REF ALLELE (ONLY NON REF HOMS CONSIDERED) & NONSENSE COUNT 5X MISSENSE
# Model 1: Treatment only
model_treatment <- coxme(Surv(extinction_gen, status) ~ Treatment + (1|rep),
                         data = survival_data)

# Model 2: Mutation load nonsense only
model_load_rec_5Weight <- coxme(Surv(extinction_gen, status) ~ load_rec_5Weight + (1|rep),
                                data = survival_data)
summary(model_load_rec_5Weight)
# Model 3: Treatment + Mutation load
model_load_rec_5Weight_both <- coxme(Surv(extinction_gen, status) ~ Treatment + load_rec_5Weight + (1|rep),
                                     data = survival_data)
summary(model_load_rec_5Weight_both)
# Extract log-likelihoods
loglik_model_treatment <- as.numeric(logLik(model_treatment))
loglik_model_load_rec_5Weight <- as.numeric(logLik(model_load_rec_5Weight))
loglik_model_load_rec_5Weight_both <- as.numeric(logLik(model_load_rec_5Weight_both))

# Number of fixed effects parameters
n_fixed_treatment <- length(fixef(model_treatment))
n_fixed_load <- length(fixef(model_load_rec_5Weight))
n_fixed_both <- length(fixef(model_load_rec_5Weight_both))

# Number of random effects variance parameters
n_random_treatment <- length(model_treatment$variance)
n_random_load <- length(model_load_rec_5Weight$variance)
n_random_both <- length(model_load_rec_5Weight_both$variance)

# Total parameters = fixed + random variance parameters
npar_model_treatment <- n_fixed_treatment + n_random_treatment
npar_model_load_rec_5Weight <- n_fixed_load + n_random_load
npar_model_load_rec_5Weight_both <- n_fixed_both + n_random_both

# Calculate AIC
aic_model_treatment <- -2 * loglik_model_treatment + 2 * npar_model_treatment
aic_model_load_rec_5Weight <- -2 * loglik_model_load_rec_5Weight + 2 * npar_model_load_rec_5Weight
aic_model_load_rec_5Weight_both <- -2 * loglik_model_load_rec_5Weight_both + 2 * npar_model_load_rec_5Weight_both

# Create comparison table
model_comparison_load_rec_5Weight <- data.frame(
  Load_scenario = "load_rec_5Weight",
  Model = c("Treatment", "Mutation Load", "Treatment + Load"),
  LogLik = c(loglik_model_treatment, loglik_model_load_rec_5Weight, loglik_model_load_rec_5Weight_both),
  n_params = c(npar_model_treatment, npar_model_load_rec_5Weight, npar_model_load_rec_5Weight_both),
  AIC = c(aic_model_treatment, aic_model_load_rec_5Weight, aic_model_load_rec_5Weight_both)
)

print(model_comparison_load_rec_5Weight)




###### 8) PARTIAL RECESSIVENESS (NON REF HOMS COUNT DOUBLE HET)  & NONSENSE COUNT 5X MISSENSE
# Model 1: Treatment only
model_treatment <- coxme(Surv(extinction_gen, status) ~ Treatment + (1|rep),
                         data = survival_data)

# Model 2: Mutation load nonsense only
model_load_parrec_5Weight <- coxme(Surv(extinction_gen, status) ~ load_parrec_5Weight + (1|rep),
                                   data = survival_data)
summary(model_load_parrec_5Weight)
# Model 3: Treatment + Mutation load
model_load_parrec_5Weight_both <- coxme(Surv(extinction_gen, status) ~ Treatment + load_parrec_5Weight + (1|rep),
                                        data = survival_data)
summary(model_load_parrec_5Weight_both)
# Extract log-likelihoods using logLik()
loglik_model_treatment <- as.numeric(logLik(model_treatment))
loglik_model_load_parrec_5Weight <- as.numeric(logLik(model_load_parrec_5Weight))
loglik_model_load_parrec_5Weight_both <- as.numeric(logLik(model_load_parrec_5Weight_both))

# Count number of fixed effects parameters
n_fixed_treatment <- length(fixef(model_treatment))
n_fixed_load <- length(fixef(model_load_parrec_5Weight))
n_fixed_both <- length(fixef(model_load_parrec_5Weight_both))

# Count number of random effects variance parameters
n_random_treatment <- length(model_treatment$variance)
n_random_load <- length(model_load_parrec_5Weight$variance)
n_random_both <- length(model_load_parrec_5Weight_both$variance)

# Total number of parameters = fixed effects + random effects variances
npar_model_treatment <- n_fixed_treatment + n_random_treatment
npar_model_load_parrec_5Weight <- n_fixed_load + n_random_load
npar_model_load_parrec_5Weight_both <- n_fixed_both + n_random_both

# Calculate AIC manually
aic_model_treatment <- -2 * loglik_model_treatment + 2 * npar_model_treatment
aic_model_load_parrec_5Weight <- -2 * loglik_model_load_parrec_5Weight + 2 * npar_model_load_parrec_5Weight
aic_model_load_parrec_5Weight_both <- -2 * loglik_model_load_parrec_5Weight_both + 2 * npar_model_load_parrec_5Weight_both

# Make a comparison table
model_comparison_load_parrec_5Weight <- data.frame(
  Load_scenario = "load_parrec_5Weight",
  Model = c("Treatment", "Mutation Load", "Treatment + Load"),
  LogLik = c(loglik_model_treatment, loglik_model_load_parrec_5Weight, loglik_model_load_parrec_5Weight_both),
  n_params = c(npar_model_treatment, npar_model_load_parrec_5Weight, npar_model_load_parrec_5Weight_both),
  AIC = c(aic_model_treatment, aic_model_load_parrec_5Weight, aic_model_load_parrec_5Weight_both)
)

print(model_comparison_load_parrec_5Weight)
      
      
      
      

###### 9) CO-DOMINANCE (NON REF HOMS COUNT DOUBLE HET) & NONSENSE COUNT 5X MISSENSE
# Model 1: Treatment only
model_treatment <- coxme(Surv(extinction_gen, status) ~ Treatment + (1|rep),
                         data = survival_data)

# Model 2: Mutation load nonsense only
model_load_codom_5Weight <- coxme(Surv(extinction_gen, status) ~ load_codom_5Weight + (1|rep),
                                  data = survival_data)
summary(model_load_codom_5Weight)
# Model 3: Treatment + Mutation load
model_load_codom_5Weight_both <- coxme(Surv(extinction_gen, status) ~ Treatment + load_codom_5Weight + (1|rep),
                                       data = survival_data)
summary(model_load_codom_5Weight_both)
# Extract log-likelihoods
loglik_model_treatment <- as.numeric(logLik(model_treatment))
loglik_model_load_codom_5Weight <- as.numeric(logLik(model_load_codom_5Weight))
loglik_model_load_codom_5Weight_both <- as.numeric(logLik(model_load_codom_5Weight_both))

# Number of fixed effects parameters
n_fixed_treatment <- length(fixef(model_treatment))
n_fixed_load <- length(fixef(model_load_codom_5Weight))
n_fixed_both <- length(fixef(model_load_codom_5Weight_both))

# Number of random effects variance parameters
n_random_treatment <- length(model_treatment$variance)
n_random_load <- length(model_load_codom_5Weight$variance)
n_random_both <- length(model_load_codom_5Weight_both$variance)

# Total parameters = fixed + random variance parameters
npar_model_treatment <- n_fixed_treatment + n_random_treatment
npar_model_load_codom_5Weight <- n_fixed_load + n_random_load
npar_model_load_codom_5Weight_both <- n_fixed_both + n_random_both

# Calculate AIC
aic_model_treatment <- -2 * loglik_model_treatment + 2 * npar_model_treatment
aic_model_load_codom_5Weight <- -2 * loglik_model_load_codom_5Weight + 2 * npar_model_load_codom_5Weight
aic_model_load_codom_5Weight_both <- -2 * loglik_model_load_codom_5Weight_both + 2 * npar_model_load_codom_5Weight_both

# Create comparison table
model_comparison_load_codom_5Weight <- data.frame(
  Load_scenario = "load_codom_5Weight",
  Model = c("Treatment", "Mutation Load", "Treatment + Load"),
  LogLik = c(loglik_model_treatment, loglik_model_load_codom_5Weight, loglik_model_load_codom_5Weight_both),
  n_params = c(npar_model_treatment, npar_model_load_codom_5Weight, npar_model_load_codom_5Weight_both),
  AIC = c(aic_model_treatment, aic_model_load_codom_5Weight, aic_model_load_codom_5Weight_both)
)

print(model_comparison_load_codom_5Weight)





###### aggregate the comparions of different models:

all_comparisons <- rbind(
  model_comparison_load_rec_noWeight,
  model_comparison_load_parrec_noWeight,
  model_comparison_load_codom_noWeight,
  model_comparison_load_rec_2Weight,
  model_comparison_load_parrec_2Weight,
  model_comparison_load_codom_2Weight,
  model_comparison_load_rec_5Weight,
  model_comparison_load_parrec_5Weight,
  model_comparison_load_codom_5Weight
)

all_comparisons <- all_comparisons %>%
  group_by(Load_scenario) %>%
  mutate(
    Delta_AIC = AIC - min(AIC),
    AIC_weight = exp(-0.5 * Delta_AIC) / sum(exp(-0.5 * Delta_AIC))
  ) %>%
  ungroup()

print(all_comparisons, n=27)
print(
  all_comparisons %>%
    mutate(
      AIC = format(AIC, digits = 4),
      LogLik = format(LogLik, digits = 4),
      Delta_AIC = format(Delta_AIC, digits = 4),
      AIC_weight = format(AIC_weight, digits = 4)
    ),
  n = 27
) %>% 
  select(-1,-2) %>% 
  print(n=27)

###########
##########
###########  PLOTS


survival_data <- survival_data %>%
  mutate(mutation_bin = ntile(mean_ratio_non_syn, 2)) %>%
  mutate(mutation_bin = factor(mutation_bin, labels = c("Low", "High")))


surv_obj <- Surv(time = survival_data$extinction_gen, event = survival_data$status)

fit_treatment <- survfit(surv_obj ~ Treatment, data = survival_data)
plot1 <- ggsurvplot(fit_treatment, data = survival_data,
                    title = "Survival by Treatment",
                    xlab = "Generation", ylab = "Survival Probability")

fit_mutload <- survfit(surv_obj ~ mutation_bin, data = survival_data)
plot2 <- ggsurvplot(fit_mutload, data = survival_data,
                    title = "Survival by Mutation Load",
                    xlab = "Generation", ylab = "Survival Probability")

ggarrange(plot1,
          plot2)









dd %>% 
  ggplot( aes( x= as.factor(gen), y=surviving, group= rep)) +
  geom_line(aes(col=Treatment)) +
  theme_bw() +
  labs(title="Lumley 2015 data - mono / poly lines",
       x= "generation of inbreeding",
       y="proportion of replicates surviving")
