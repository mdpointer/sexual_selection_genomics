#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 10G
##SBATCH --cpus-per-task 4
#SBATCH -t 150:00:00                               # Set time limit 
#SBATCH --job-name=ss.ARRAY_get_alelle_freqs_per_pop_per_snp_cat_intergenicRegions       #Job name
#SBATCH -o ss.ARRAY_get_alelle_freqs_per_pop_per_snp_cat_intergenicRegions.out           #Standard output log
#SBATCH -e ss.ARRAY_get_alelle_freqs_per_pop_per_snp_cat_intergenicRegions.err           #Standard error log

#SBATCH --array=0-22

# # set up the environment
module add SIFT4G/2.0.0
module add bcftools
module add htslib/1.15.1

MASTER=~/ss

INTERGENIC_SNP_VCF=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_GA1_mono_poly__intergenicRegions.vcf.gz


OUTPUT_FOLDER=${MASTER}/allele_Fs_intergenicRegions

PATH_TO_POPFILES=${MASTER}/population_files/182sample_popfiles
POPS=(GA1 monoall polyall GA1_males GA1_females monoA monoA_males monoA_females monoB monoB_males monoB_females monoC monoC_males monoC_females polyA polyA_males polyA_females polyB polyB_males polyB_females polyC polyC_males polyC_females)

POP_FILE=${PATH_TO_POPFILES}/ss.samples_${POPS[${SLURM_ARRAY_TASK_ID}]}.txt





# ## FOR EACH VARIANT CATEGORY, GET THE ALT ALLELE FREQUENCY, ADD SOME IDENTIFYING COLUMNS, INCLUDING A CHROM_POS COLUMN, THEN SORT ON THAT COLUMN

# INTERGENIC
bcftools view -S ${POP_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${INTERGENIC_SNP_VCF} | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | awk -v pop="${POPS[${SLURM_ARRAY_TASK_ID}]}" -v VCF_info="intergenic" -v chrom="autosomes" 'BEGIN { OFS="\t" } { combined_pos = $1 "_" $2; print $0, pop, VCF_info, chrom, combined_pos }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_INTERGENIC_autosomes_alleleFs.tsv
sort -k9,9 ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_INTERGENIC_autosomes_alleleFs.tsv > ${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_INTERGENIC_autosomes_alleleFs.tsv

bcftools view -S ${POP_FILE} -r LGX ${INTERGENIC_SNP_VCF} | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | awk -v pop="${POPS[${SLURM_ARRAY_TASK_ID}]}" -v VCF_info="intergenic" -v chrom="X" 'BEGIN { OFS="\t" } { combined_pos = $1 "_" $2; print $0, pop, VCF_info, chrom, combined_pos }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_INTERGENIC_X_alleleFs.tsv
sort -k9,9 ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_INTERGENIC_X_alleleFs.tsv > ${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_INTERGENIC_X_alleleFs.tsv


# remove unsorted files to save space
rm -f ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_INTERGENIC_autosomes_alleleFs.tsv ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_INTERGENIC_X_alleleFs.tsv 






###  FOR EACH FILE CREATED ABOVE, JOIN IT WITH A FILE CONTAINING THE MINOR ALLELE IN GA1 AT EVERY LOCUS

# # INTERGENIC
join -1 9 -2 6 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.6 2.5 -e "NA" -a 1 -t $'\t' \
${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_INTERGENIC_autosomes_alleleFs.tsv \
~/ss/allele_Fs_intergenicRegions/sorted_GA1_minor_alleles_INTERGENIC_autosomes.txt | \
awk 'BEGIN { OFS="\t" } { if ($4 == $11) $12 = $5; else $12 = 1 - $5; print }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_INTERGENIC_autosomes_minor_alleleFs.tsv

join -1 9 -2 6 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.6 2.5 -e "NA" -a 1 -t $'\t' \
${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_INTERGENIC_X_alleleFs.tsv \
~/ss/allele_Fs_intergenicRegions/sorted_GA1_minor_alleles_INTERGENIC_X.txt | \
awk 'BEGIN { OFS="\t" } { if ($4 == $11) $12 = $5; else $12 = 1 - $5; print }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_INTERGENIC_X_minor_alleleFs.tsv





echo "End of script"

