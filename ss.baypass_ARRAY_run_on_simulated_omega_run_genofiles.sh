#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH -t 48:00:00                               # Set time limit to 36 hours
#SBATCH --job-name=ss.baypass_ARRAY_run_on_simulated_omega_run_genofiles        #Job name
#SBATCH --mem 20G
#SBATCH --cpus-per-task 8
#SBATCH -o ss.baypass_ARRAY_run_on_simulated_omega_run_genofiles_%A_%a.out                    #Standard output log
#SBATCH -e ss.baypass_ARRAY_run_on_simulated_omega_run_genofiles_%A_%a.err                     #Standard error log

#SBATCH --array=0

# Set up environment
module load baypass/2.4
echo "baypass loaded"

MASTER=~/ss


#### to run baypass on the simulated genofile you need:
# simulated genofile
# number of pops
# covariate file if you're using one
# output prefix

NPOP=(6 6)

SIM_GENOFILE_PATH=~/ss/baypass/simulated_genofiles_omega_runs
SIM_GENOFILES=(
    "${SIM_GENOFILE_PATH}/G.MbiasedFbiased_SIM"
    "${SIM_GENOFILE_PATH}/G.MonoPoly_SIM"
)

OMEGAFILE_PATH=~/ss/baypass/omega_matrices
OMEGAFILES=(
    "${OMEGAFILE_PATH}/omega_run_6pop_nonEXONIC10kb_onlyMbFb_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt_mat_omega.out"
    "${OMEGAFILE_PATH}/omega_run_6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt_mat_omega.out"
)

CONTRASTFILE_PATH=~/ss/data/baypass_contrast_files
CONTRASTFILES=(
    "${CONTRASTFILE_PATH}/ss.baypass_MbFb_contrast.txt"
    "${CONTRASTFILE_PATH}/ss.baypass_MonoPoly_contrast.txt"
)

OUTPREFIX_PATH=~/ss/baypass/POD/simulated_genofiles_data_runs
OUTPREFIX=(
    "${OUTPREFIX_PATH}/ss.baypass_MbiasedFbiased_SIM"
    "${OUTPREFIX_PATH}/ss.baypass_MonoPoly_SIM"
)

echo "All input files loaded - about to call Baypass"

# run Baypass independently on all available SIMULATED genetype matrices

baypass \
-npop ${NPOP[${SLURM_ARRAY_TASK_ID}]} \
-gfile ${SIM_GENOFILES[${SLURM_ARRAY_TASK_ID}]} \
-omegafile ${OMEGAFILES[${SLURM_ARRAY_TASK_ID}]} \
-efile ${CONTRASTFILES[${SLURM_ARRAY_TASK_ID}]} \
-contrastfile ${CONTRASTFILES[${SLURM_ARRAY_TASK_ID}]} \
-auxmodel \
-outprefix ${OUTPREFIX[${SLURM_ARRAY_TASK_ID}]} \
-nthreads 8 \
-seed 1000