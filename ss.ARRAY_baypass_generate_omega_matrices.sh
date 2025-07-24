#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 20G
#SBATCH --cpus-per-task 8
#SBATCH -t 100:00:00                               # Set time limit to 100 hours
#SBATCH --job-name=ss.ARRAY_baypass_generate_omega_matrices        #Job name
#SBATCH -o ss.ARRAY_baypass_generate_omega_matrices_%A_%a.out     #Standard output log
#SBATCH -e ss.ARRAY_baypass_generate_omega_matrices_%A_%a.err     #Standard error log

#SBATCH --array=0

# Set up environment
module load baypass/2.4

MASTER=~/ss

# Define an array with filenames
GENOFILES=(
    "6pop_nonEXONIC10kb_onlyMbFb_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt"
    "6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt"
)

echo "GENOFILE FOR THIS RUN IS: ${GENOFILES[${SLURM_ARRAY_TASK_ID}]}"



# set up for Baypass
NPOP=6
mkdir -p ${MASTER}/baypass/omega_matrices

# run Baypass independently on all available genetype matrices
baypass \
-npop ${NPOP} \
-gfile ${MASTER}/data/genofiles_for_omega_run/${GENOFILES[${SLURM_ARRAY_TASK_ID}]} \
-outprefix ${MASTER}/baypass/omega_matrices/omega_run_${GENOFILES[${SLURM_ARRAY_TASK_ID}]} \
-nthreads 8

echo "Baypass run to generate omega matrix from	${GENOFILES[${SLURM_ARRAY_TASK_ID}]} ... COMPLETED"