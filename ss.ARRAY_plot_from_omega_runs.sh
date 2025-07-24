#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH -t 48:00:00                               # Set time limit to 36 hours
#SBATCH --job-name=ARRAY_plot_from_omega_runs        #Job name
#SBATCH --mem 20G
#SBATCH -o ARRAY_plot_from_omega_runs_%A_%a.out                    #Standard output log
#SBATCH -e ARRAY_plot_from_omega_runs_%A_%a.err                     #Standard error log

#SBATCH --array=0-7


# set u[ the environment
module load R/4.3.1

NPOP=(32 32 32 32 2 2 2 2)


# load the XtX output of omega runs
OMEGA_XTX=(
  "omega_run_32pop_allREG_onlyHL_allSCAFF_MAF003_MAXMISS2_LD0350_genotypes.txt_summary_pi_xtx.out"
  "omega_run_32pop_allREG_onlyHL_onlyCHROM_MAF003_MAXMISS2_LD0350_genotypes.txt_summary_pi_xtx.out"
  "omega_run_32pop_nonEXONIC10kb_onlyHL_allSCAFF_MAF003_MAXMISS2_LD0350_genotypes.txt_summary_pi_xtx.out" 
  "omega_run_32pop_nonEXONIC10kb_onlyHL_onlyCHROM_MAF003_MAXMISS2_LD0350_genotypes.txt_summary_pi_xtx.out"
  "omega_run_2pop_allREG_onlyHL_allSCAFF_MAF003_MAXMISS2_LD0350_genotypes.txt_summary_pi_xtx.out"
  "omega_run_2pop_allREG_onlyHL_onlyCHROM_MAF003_MAXMISS2_LD0350_genotypes.txt_summary_pi_xtx.out"
  "omega_run_2pop_nonEXONIC10kb_onlyHL_allSCAFF_MAF003_MAXMISS2_LD0350_genotypes.txt_summary_pi_xtx.out"
  "omega_run_2pop_nonEXONIC10kb_onlyHL_onlyCHROM_MAF003_MAXMISS2_LD0350_genotypes.txt_summary_pi_xtx.out"
)


# load SNP position files for each run
SNP_POS=(
  "allREG_onlyHL_allSCAFF_MAF003_MAXMISS2_LD0350.txt"
  "allREG_onlyHL_onlyCHROM_MAF003_MAXMISS2_LD0350.txt"
  "nonEXONIC10kb_onlyHL_allSCAFF_MAF003_MAXMISS2_LD0350.txt"
  "nonEXONIC10kb_onlyHL_onlyCHROM_MAF003_MAXMISS2_LD0350.txt"
  "allREG_onlyHL_allSCAFF_MAF003_MAXMISS2_LD0350.txt"
  "allREG_onlyHL_onlyCHROM_MAF003_MAXMISS2_LD0350.txt"
  "nonEXONIC10kb_onlyHL_allSCAFF_MAF003_MAXMISS2_LD0350.txt"
  "nonEXONIC10kb_onlyHL_onlyCHROM_MAF003_MAXMISS2_LD0350.txt"
)


# load the XtX outputs of runs using simulated genofiles from omega run genofiles (to use for thresholds
POD_XTX=(
    "32pop_allREG_allSCAFF_POD_summary_pi_xtx.out"
    "32pop_allREG_onlyCHROM_POD_summary_pi_xtx.out"
    "32pop_nonEXONIC10kb_allSCAFF_POD_summary_pi_xtx.out"
    "32pop_nonEXONIC10kb_onlyCHROM_POD_summary_pi_xtx.out"
    "2pop_allREG_allSCAFF_POD_summary_pi_xtx.out"
    "2pop_allREG_onlyCHROM_POD_summary_pi_xtx.out"
    "2pop_nonEXONIC10kb_allSCAFF_POD_summary_pi_xtx.out"
    "2pop_nonEXONIC10kb_onlyCHROM_POD_summary_pi_xtx.out"
)



# arguments to pass to the R script
arg1=${NPOP[${SLURM_ARRAY_TASK_ID}]}
arg2=${OMEGA_XTX[${SLURM_ARRAY_TASK_ID}]}
arg3=${SNP_POS[${SLURM_ARRAY_TASK_ID}]}
arg4=${POD_XTX[${SLURM_ARRAY_TASK_ID}]}


echo " ${arg1} ${arg2} ${arg3} ${arg4} "

Rscript ~/dispersal/scripts/R_baypass/ARRAY_plot_from_omega_runs.R "${arg1}" "${arg2}" "${arg3}" "${arg4}"