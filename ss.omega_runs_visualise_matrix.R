library(tidyverse)
library(stats)
library(corrplot)
library(qvalue)
#rm(list=ls())

######### Plot omega matrix for Male-biased / Female-biased

### visualise covariance matrix, convert covariance to correlation with cov2cor()
omegafile_filename <- "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/baypass/omega_matrices/omega_run_6pop_nonEXONIC10kb_onlyMbFb_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt_mat_omega.out"
covmat <- read_delim( omegafile_filename, delim="    ", col_names = F) %>%  # I think the delim is \t, if it doesn't work then try others
  select(c(paste(rep("X", 6), seq(2, 12, 2), sep=""))) %>% 
  rename(M_biased_A=X2,
         M_biased_B=X4,
         M_biased_C=X6,
         F_biased_A=X8,
         F_biased_B=X10,
         F_biased_C=X12)
covmat <- covmat %>% as.matrix()
cormat <- stats::cov2cor(covmat)
rownames(cormat) <- c("M_biased_A", "M_biased_B", "M_biased_C", "F_biased_A", "F_biased_B", "F_biased_C")
corrplot(cormat, method='color', type='lower', diag=T)



######### Plot omega matrix for Mono / Poly

### visualise covariance matrix, convert covariance to correlation with cov2cor()
omegafile_filename <- "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/baypass/omega_matrices/omega_run_6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt_mat_omega.out"
covmat <- read_delim( omegafile_filename, delim="    ", col_names = F) %>%  # I think the delim is \t, if it doesn't work then try others
  select(c(paste(rep("X", 6), seq(2, 12, 2), sep=""))) %>% 
  rename(mono_A=X2,
         mono_B=X4,
         mono_C=X6,
         poly_A=X8,
         poly_B=X10,
         poly_C=X12)
covmat <- covmat %>% as.matrix()
cormat <- stats::cov2cor(covmat)
rownames(cormat) <- c("mono_A", "mono_B", "mono_C", "poly_A", "poly_B", "poly_C")
corrplot(cormat, method='color', type='lower', diag=T)






