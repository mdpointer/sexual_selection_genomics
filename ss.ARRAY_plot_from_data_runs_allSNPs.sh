#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-24-96                          #Which queue to use
#SBATCH -t 48:00:00                               # Set time limit to 36 hours
#SBATCH --job-name=ss.ARRAY_plot_from_data_runs_allSNPs        #Job name
#SBATCH --mem 20G
#SBATCH -o ss.ARRAY_plot_from_data_runs_allSNPs_%A_%a.out                    #Standard output log
#SBATCH -e ss.ARRAY_plot_from_data_runs_allSNPs_%A_%a.err                     #Standard error log

#SBATCH --array=0-1


# set up the environment
module load R/4.3.1

NPOP=(6 6)



##
## specify arrays of output file names NOTE - file paths are provided in the Rscript called below

# load the XtX output of data runs
XTX=(
  ss_baypass_datarun_6pop_MbFb_onlyAUTOSOMES_allREG_allSNPs_summary_pi_xtx.out
  ss_baypass_datarun_6pop_MonoPoly_onlyAUTOSOMES_allREG_allSNPs_summary_pi_xtx.out
)

# load the BF output of data runs
BF=(
  ss_baypass_datarun_6pop_MbFb_onlyAUTOSOMES_allREG_allSNPs_summary_betai.out
  ss_baypass_datarun_6pop_MonoPoly_onlyAUTOSOMES_allREG_allSNPs_summary_betai.out
)

# load the C2 output of data runs
C2=(
  ss_baypass_datarun_6pop_MbFb_onlyAUTOSOMES_allREG_allSNPs_summary_contrast.out
  ss_baypass_datarun_6pop_MonoPoly_onlyAUTOSOMES_allREG_allSNPs_summary_contrast.out
)

# load SNP position files for each run
SNP_POS=(
  allREGION_onlyMbFb_onlyLGautosomes_MAF003_noLDprune.txt
  allREGION_onlyMonoPoly_onlyLGautosomes_MAF003_noLDprune.txt
)

# load POD run simulated xtx files
POD_XTX=(
  ss.baypass_MbiasedFbiased_SIM_summary_pi_xtx.out
  ss.baypass_MonoPoly_SIM_summary_pi_xtx.out
)

# load POD run simulated C2 files
POD_C2=(
  ss.baypass_MbiasedFbiased_SIM_summary_contrast.out
  ss.baypass_MonoPoly_SIM_summary_contrast.out
)


# arguments to pass to the R script
arg1=${NPOP[${SLURM_ARRAY_TASK_ID}]}
arg2=${XTX[${SLURM_ARRAY_TASK_ID}]}
arg3=${BF[${SLURM_ARRAY_TASK_ID}]}
arg4=${C2[${SLURM_ARRAY_TASK_ID}]}
arg5=${SNP_POS[${SLURM_ARRAY_TASK_ID}]}
arg6=${POD_XTX[${SLURM_ARRAY_TASK_ID}]}
arg7=${POD_C2[${SLURM_ARRAY_TASK_ID}]}

echo " ${arg1} ${arg2} ${arg3} ${arg4} ${arg5} ${arg6} ${arg7}"

Rscript ~/ss/scripts/baypass/ss.ARRAY_plot_from_data_runs_allSNPs.R "${arg1}" "${arg2}" "${arg3}" "${arg4}" "${arg5}" "${arg6}" "${arg7}"