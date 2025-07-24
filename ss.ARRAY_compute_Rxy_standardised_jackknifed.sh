#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 5G
##SBATCH --cpus-per-task 4
#SBATCH -t 150:00:00                               # Set time limit 
#SBATCH --job-name=ss.ARRAY_compute_Rxy_standardised_jackknifed       #Job name
#SBATCH -o ss.ARRAY_compute_Rxy_standardised_jackknifed.out           #Standard output log
#SBATCH -e ss.ARRAY_compute_Rxy_standardised_jackknifed.err           #Standard error log

#SBATCH --array=0-23


# # set up the environment
module add R/4.3.1


POPX=(
    polyall
    polyall
    polyall
    polyall
    polyall
    polyall
    polyall
    polyall
    monoall
    monoall
    monoall
    monoall
    polyall
    polyall
    polyall
    polyall
    polyall
    polyall
    polyall
    polyall
    monoall
    monoall
    monoall
    monoall
)

POPX_codingRegions_filename=(
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_SYNONYMOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_MISSENSE_TOLERATED_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_MISSENSE_DELETERIOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_NONSENSE_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_SYNONYMOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_MISSENSE_TOLERATED_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_MISSENSE_DELETERIOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_NONSENSE_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_SYNONYMOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_MISSENSE_TOLERATED_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_MISSENSE_DELETERIOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_NONSENSE_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_SYNONYMOUS_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_MISSENSE_TOLERATED_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_MISSENSE_DELETERIOUS_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_NONSENSE_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_SYNONYMOUS_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_MISSENSE_TOLERATED_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_MISSENSE_DELETERIOUS_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/polyall_NONSENSE_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_SYNONYMOUS_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_MISSENSE_TOLERATED_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_MISSENSE_DELETERIOUS_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_NONSENSE_X_minor_alleleFs.tsv
)

POPX_intergenicRegions_filename=(
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/polyall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_X_minor_alleleFs.tsv
)

POPY=(
    monoall
    monoall
    monoall
    monoall
    GA1
    GA1
    GA1
    GA1
    GA1
    GA1
    GA1
    GA1
    monoall
    monoall
    monoall
    monoall
    GA1
    GA1
    GA1
    GA1
    GA1
    GA1
    GA1
    GA1
)

POPY_codingRegions_filename=(
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_SYNONYMOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_MISSENSE_TOLERATED_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_MISSENSE_DELETERIOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_NONSENSE_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_SYNONYMOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_MISSENSE_TOLERATED_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_MISSENSE_DELETERIOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_NONSENSE_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_SYNONYMOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_MISSENSE_TOLERATED_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_MISSENSE_DELETERIOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_NONSENSE_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_SYNONYMOUS_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_MISSENSE_TOLERATED_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_MISSENSE_DELETERIOUS_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/monoall_NONSENSE_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_SYNONYMOUS_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_MISSENSE_TOLERATED_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_MISSENSE_DELETERIOUS_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_NONSENSE_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_SYNONYMOUS_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_MISSENSE_TOLERATED_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_MISSENSE_DELETERIOUS_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_codingRegions_SIFToutput/GA1_NONSENSE_X_minor_alleleFs.tsv
)


POPY_intergenicRegions_filename=(
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_autosomes_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/monoall_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_X_minor_alleleFs.tsv
    ~/ss/allele_Fs_intergenicRegions/GA1_INTERGENIC_X_minor_alleleFs.tsv

)

VARIANT_CAT=(
    synonymous
    missense_tolerated
    missense_deleterious
    nonsense
    synonymous
    missense_tolerated
    missense_deleterious
    nonsense
    synonymous
    missense_tolerated
    missense_deleterious
    nonsense
    synonymous
    missense_tolerated
    missense_deleterious
    nonsense
    synonymous
    missense_tolerated
    missense_deleterious
    nonsense
    synonymous
    missense_tolerated
    missense_deleterious
    nonsense
)

CHROM=(
    autosomes
    autosomes
    autosomes
    autosomes
    autosomes
    autosomes
    autosomes
    autosomes
    autosomes
    autosomes
    autosomes
    autosomes
    X
    X
    X
    X
    X
    X
    X
    X
    X
    X
    X
    X
)

arg1=${POPX[${SLURM_ARRAY_TASK_ID}]}
arg2=${POPX_codingRegions_filename[${SLURM_ARRAY_TASK_ID}]}   # population X allele frequency file name
arg3=${POPX_intergenicRegions_filename[${SLURM_ARRAY_TASK_ID}]}
arg4=${POPY[${SLURM_ARRAY_TASK_ID}]}
arg5=${POPY_codingRegions_filename[${SLURM_ARRAY_TASK_ID}]}   # population Y allele frequency file name
arg6=${POPY_intergenicRegions_filename[${SLURM_ARRAY_TASK_ID}]}
arg7=${VARIANT_CAT[${SLURM_ARRAY_TASK_ID}]} 
arg8=${CHROM[${SLURM_ARRAY_TASK_ID}]}     # autosomes / X

Rscript ~/ss/scripts/Rxy/ss_compute_Rxy_standardised_jackknifed.R "${arg1}" "${arg2}" "${arg3}" "${arg4}" "${arg5}" "${arg6}" "${arg7}" "${arg8}"

