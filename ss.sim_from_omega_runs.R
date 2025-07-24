library(tidyverse)
#install.packages("mvtnorm")
library(mvtnorm)

setwd("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/baypass/omega_runs/sim_genofiles_baypass/sim_output")
# to do the sim for 1 run you need:
# omega matrix
# pi.beta.coefs - if using a covariable
# genotype matrix


# run the Baypass utils script to get functions for simulating POD
source("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/baypass/omega_runs/baypass_utils.R")


# vector of empirical omega matrix files
emp_omegafiles <- c(
  "omega_run_6pop_nonEXONIC10kb_onlyMbFb_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt_mat_omega.out",
  "omega_run_6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt_mat_omega.out"
)
emp_omegafiles_location <- "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/baypass/omega_runs/sim_genofiles_baypass/omega_run_emp_omegafiles/"

# use the names of the empirical genofiles
emp_genofiles <- c(
 "6pop_nonEXONIC10kb_onlyMbFb_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt",
 "6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt"
)
emp_genofiles_location <- "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/baypass/omega_runs/sim_genofiles_baypass/omega_run_emp_genofiles/"


runs <- c("MbiasedFbiased", "MonoPoly")

# loop over different runs
for (i in 1:length(emp_genofiles)) {
  
  # define and load the omegafile for the focal run
  omega_filename <- paste(emp_omegafiles_location, emp_omegafiles[i], sep="")
  emp_omega_matrix <- as.matrix(read.table(omega_filename))
  
  # define and load the genofile for the focal run
  genofile_filename <- paste(emp_genofiles_location, emp_genofiles[i], sep="")
  # emp_genofile <- 
  
  # pi.beta.coef = read.table("baypass_output_191B1K+53KS_in4pops_summary_beta_params.out",h=T)$Mean # required to sim covar model, but not core
  
  # simulate a neutral genotype matrix of nsnp SNPs and sample.size indivs per pop
  bp_sim <- simulate.baypass(
    omega.mat=emp_omega_matrix,
    nsnp=100000,
    beta.coef=NA,
    #beta.pi=c(1,1),
    #pop.trait= env_cov_vec,
    sample.size=24, # you can get this by running geno2yn on genofile and taking the output$NN, but for me it's just always 6 samples per pop = 12 alleles
    pi.maf=0.04,
    suffix= paste(runs[i],"_SIM", sep=""),
    remove.fixed.loci=T,
    coverage=NA
    )
  
}
  