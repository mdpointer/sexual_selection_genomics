#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 4G
##SBATCH --cpus-per-task 4
#SBATCH -t 150:00:00                               # Set time limit 
#SBATCH --job-name=ss.ALLELE_Fs_autosomes_import_and_plot       #Job name
#SBATCH -o ss.ALLELE_Fs_autosomes_import_and_plot.out           #Standard output log
#SBATCH -e ss.ALLELE_Fs_autosomes_import_and_plot.err           #Standard error log


#set up environment
module load R/4.1.2


#run the R script
Rscript ~/ss/scripts/allele_Fs/ss.ALLELE_Fs_autosomes_import_and_plot