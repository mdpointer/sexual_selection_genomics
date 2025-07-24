#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 5G
##SBATCH --cpus-per-task 4
#SBATCH -t 150:00:00                               # Set time limit 
#SBATCH --job-name=ss.ARRAY_cat_outputs_and_summarise_across_jackknife       #Job name
#SBATCH -o ss.ARRAY_cat_outputs_and_summarise_across_jackknife.out           #Standard output log
#SBATCH -e ss.ARRAY_cat_outputs_and_summarise_across_jackknife.err           #Standard error log


module add R/4.3.1

cat ~/ss/Rxy/*.tsv > ~/ss/Rxy/processed_output/Rxy_combined_jackknifes_across_comparisons.tsv

COMBINED_RXY_OUTPUT_FILE=~/ss/Rxy/processed_output/Rxy_combined_jackknifes_across_comparisons.tsv




# provide the file to the R script
arg1=${COMBINED_RXY_OUTPUT_FILE}

Rscript ~/ss/scripts/Rxy/ss_Rxy_summarise_across_bootstraps.R "${arg1}"