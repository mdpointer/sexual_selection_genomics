#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 6G
#SBATCH --cpus-per-task 8
#SBATCH -t 150:00:00                               # Set time limit 
#SBATCH --job-name=ss.baypass_onlyLGX_bothSEX_data_runs_contrast        #Job name
#SBATCH -o ss.baypass_onlyLGX_bothSEX_data_runs_contrast_%A_%a.out     #Standard output log
#SBATCH -e ss.baypass_onlyLGX_bothSEX_data_runs_contrast_%A_%a.err     #Standard error log

#SBATCH --array=0-1



module load baypass/2.4

mkdir -p ~/ss/baypass/data_runs_C2_BF/combined_sex_X



# Npop
NPOP=(6 6)



# genotype matrices
GENOFILES=(
    ~/ss/data/genofiles_for_data_run/6pop_allREGION_onlyMbFb_bothSEX_onlyLGX_MAF003_noLDprune_genotypes.txt
    ~/ss/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_genotypes.txt
)



# Previously generated Omega matrix files. Made using various filterings of full SNPset
OMEGAFILES=(
    ~/ss/baypass/omega_matrices/combined_sex_X/omega_run_6pop_nonEXONIC10kb_onlyMbFb_bothSEX_onlyX_MAF003_MAXMISS1_LD0350_genotypes.txt_mat_omega.out
    ~/ss/baypass/omega_matrices/combined_sex_X/omega_run_6pop_nonEXONIC10kb_onlyMonoPoly_bothSEX_onlyX_MAF003_MAXMISS1_LD0350_genotypes.txt_mat_omega.out
)


# binary population covariates
CONTRASTFILES=(
    ~/ss/data/baypass_contrast_files/ss.baypass_MbFb_contrast.txt
    ~/ss/data/baypass_contrast_files/ss.baypass_MonoPoly_contrast.txt
)



# OUPREFIX
OUTPREFIXES=(
    ~/ss/baypass/data_runs_C2_BF/combined_sex_X/ss_baypass_datarun_6pop_MbFb_bothSEX_onlyX_allREG
    ~/ss/baypass/data_runs_C2_BF/combined_sex_X/ss_baypass_datarun_6pop_MonoPoly_bothSEX_onlyX_allREG
)





##
## RUN BAYPASS
##

echo "Beginning run ${OUTPREFIXES[${SLURM_ARRAY_TASK_ID}]} "

baypass \
-npop ${NPOP[${SLURM_ARRAY_TASK_ID}]} \
-gfile ${GENOFILES[${SLURM_ARRAY_TASK_ID}]} \
-omegafile ${OMEGAFILES[${SLURM_ARRAY_TASK_ID}]} \
-efile ${CONTRASTFILES[${SLURM_ARRAY_TASK_ID}]} \
-contrastfile ${CONTRASTFILES[${SLURM_ARRAY_TASK_ID}]} \
-auxmodel \
-outprefix ${OUTPREFIXES[${SLURM_ARRAY_TASK_ID}]} \
-nthreads 8 \
-seed 1000

echo "Finished run ${OUTPREFIXES[${SLURM_ARRAY_TASK_ID}]} "