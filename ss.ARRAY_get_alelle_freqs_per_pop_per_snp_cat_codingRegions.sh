#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 10G
##SBATCH --cpus-per-task 4
#SBATCH -t 150:00:00                               # Set time limit 
#SBATCH --job-name=ss.ARRAY_get_alelle_freqs_per_pop_per_snp_cat_codingRegions       #Job name
#SBATCH -o ss.ARRAY_get_alelle_freqs_per_pop_per_snp_cat_codingRegions.out           #Standard output log
#SBATCH -e ss.ARRAY_get_alelle_freqs_per_pop_per_snp_cat_codingRegions.err           #Standard error log

#SBATCH --array=0-22

# # set up the environment
module add SIFT4G/2.0.0
module add bcftools
module add htslib/1.15.1

MASTER=~/ss

MISSENSE_DELETERIOUS_VCF_HIGH_CONFIDENCE=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_SIFT_MISSENSE_DELETERIOUS_high_confidence.vcf.gz
MISSENSE_TOLERATED_VCF=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_SIFT_MISSENSE_TOLERATED.vcf.gz
NONSENSE_SNP_VCF=${MASTER}/snpEff/codingRegions_snpEff_NONSENSE.vcf.gz
SYNONYMOUS_SNP_VCF=${MASTER}/snpEff/codingRegions_snpEff_SYNONYMOUS.vcf.gz

OUTPUT_FOLDER=${MASTER}/allele_Fs_codingRegions_SIFToutput

PATH_TO_POPFILES=${MASTER}/population_files/182sample_popfiles
POPS=(GA1 monoall polyall GA1_males GA1_females monoA monoA_males monoA_females monoB monoB_males monoB_females monoC monoC_males monoC_females polyA polyA_males polyA_females polyB polyB_males polyB_females polyC polyC_males polyC_females)

POP_FILE=${PATH_TO_POPFILES}/ss.samples_${POPS[${SLURM_ARRAY_TASK_ID}]}.txt





# ## FOR EACH VARIANT CATEGORY, GET THE ALT ALLELE FREQUENCY, ADD SOME IDENTIFYING COLUMNS, INCLUDING A CHROM_POS COLUMN, THEN SORT ON THAT COLUMN

# MISSENSE DEL
bcftools view -S ${POP_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${MISSENSE_DELETERIOUS_VCF_HIGH_CONFIDENCE} | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | awk -v pop="${POPS[${SLURM_ARRAY_TASK_ID}]}" -v VCF_info="missense_del" -v chrom="autosomes" 'BEGIN { OFS="\t" } { combined_pos = $1 "_" $2; print $0, pop, VCF_info, chrom, combined_pos }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_DELETERIOUS_autosomes_alleleFs.tsv
sort -k9,9 ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_DELETERIOUS_autosomes_alleleFs.tsv > ${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_DELETERIOUS_autosomes_alleleFs.tsv
bcftools view -S ${POP_FILE} -r LGX ${MISSENSE_DELETERIOUS_VCF_HIGH_CONFIDENCE} | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | awk -v pop="${POPS[${SLURM_ARRAY_TASK_ID}]}" -v VCF_info="missense_del" -v chrom="X" 'BEGIN { OFS="\t" } { combined_pos = $1 "_" $2; print $0, pop, VCF_info, chrom, combined_pos }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_DELETERIOUS_X_alleleFs.tsv
sort -k9,9 ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_DELETERIOUS_X_alleleFs.tsv > ${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_DELETERIOUS_X_alleleFs.tsv

# MISSENSE TOL
bcftools view -S ${POP_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${MISSENSE_TOLERATED_VCF} | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | awk -v pop="${POPS[${SLURM_ARRAY_TASK_ID}]}" -v VCF_info="missense_tol" -v chrom="autosomes" 'BEGIN { OFS="\t" } { combined_pos = $1 "_" $2; print $0, pop, VCF_info, chrom, combined_pos }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_TOLERATED_autosomes_alleleFs.tsv
sort -k9,9 ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_TOLERATED_autosomes_alleleFs.tsv > ${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_TOLERATED_autosomes_alleleFs.tsv
bcftools view -S ${POP_FILE} -r LGX ${MISSENSE_TOLERATED_VCF} | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | awk -v pop="${POPS[${SLURM_ARRAY_TASK_ID}]}" -v VCF_info="missense_tol" -v chrom="X" 'BEGIN { OFS="\t" } { combined_pos = $1 "_" $2; print $0, pop, VCF_info, chrom, combined_pos }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_TOLERATED_X_alleleFs.tsv
sort -k9,9 ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_TOLERATED_X_alleleFs.tsv > ${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_TOLERATED_X_alleleFs.tsv

# NONSENSE
bcftools view -S ${POP_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${NONSENSE_SNP_VCF} | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | awk -v pop="${POPS[${SLURM_ARRAY_TASK_ID}]}" -v VCF_info="nonsense" -v chrom="autosomes" 'BEGIN { OFS="\t" } { combined_pos = $1 "_" $2; print $0, pop, VCF_info, chrom, combined_pos }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_NONSENSE_autosomes_alleleFs.tsv
sort -k9,9 ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_NONSENSE_autosomes_alleleFs.tsv > ${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_NONSENSE_autosomes_alleleFs.tsv
bcftools view -S ${POP_FILE} -r LGX ${NONSENSE_SNP_VCF} | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | awk -v pop="${POPS[${SLURM_ARRAY_TASK_ID}]}" -v VCF_info="nonsense" -v chrom="X" 'BEGIN { OFS="\t" } { print $0, pop, VCF_info, chrom }'> ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_NONSENSE_X_alleleFs.tsv
sort -k9,9 ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_NONSENSE_X_alleleFs.tsv > ${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_NONSENSE_X_alleleFs.tsv

# SYNONYMOUS
bcftools view -S ${POP_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${SYNONYMOUS_SNP_VCF} | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | awk -v pop="${POPS[${SLURM_ARRAY_TASK_ID}]}" -v VCF_info="synonymous" -v chrom="autosomes" 'BEGIN { OFS="\t" } { combined_pos = $1 "_" $2; print $0, pop, VCF_info, chrom, combined_pos }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_SYNONYMOUS_autosomes_alleleFs.tsv
sort -k9,9 ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_SYNONYMOUS_autosomes_alleleFs.tsv > ${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_SYNONYMOUS_autosomes_alleleFs.tsv
bcftools view -S ${POP_FILE} -r LGX ${SYNONYMOUS_SNP_VCF} | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | awk -v pop="${POPS[${SLURM_ARRAY_TASK_ID}]}" -v VCF_info="synonymous" -v chrom="X" 'BEGIN { OFS="\t" } { combined_pos = $1 "_" $2; print $0, pop, VCF_info, chrom, combined_pos }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_SYNONYMOUS_X_alleleFs.tsv
sort -k9,9 ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_SYNONYMOUS_X_alleleFs.tsv > ${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_SYNONYMOUS_X_alleleFs.tsv


# remove unsorted files to save space
rm -f ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_DELETERIOUS_autosomes_alleleFs.tsv ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_DELETERIOUS_X_alleleFs.tsv ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_TOLERATED_autosomes_alleleFs.tsv ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_TOLERATED_X_alleleFs.tsv ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_NONSENSE_autosomes_alleleFs.tsv ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_NONSENSE_X_alleleFs.tsv ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_SYNONYMOUS_autosomes_alleleFs.tsv ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_SYNONYMOUS_autosomes_alleleFs.tsv






###  FOR EACH FILE CREATED ABOVE, JOIN IT WITH A FILE CONTAINING THE MINOR ALLELE IN GA1 AT EVERY LOCUS

# # MISSENSE DEL
join -1 9 -2 6 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.6 2.5 -e "NA" -a 1 -t $'\t' \
${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_DELETERIOUS_autosomes_alleleFs.tsv \
~/ss/allele_Fs_codingRegions_SIFToutput/sorted_GA1_minor_alleles_MISSENSE_DELETERIOUS_HIGH_CONFIDENCE_autosomes.txt | \
awk 'BEGIN { OFS="\t" } { if ($4 == $11) $12 = $5; else $12 = 1 - $5; print }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_DELETERIOUS_autosomes_minor_alleleFs.tsv

join -1 9 -2 6 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.6 2.5 -e "NA" -a 1 -t $'\t' \
${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_DELETERIOUS_X_alleleFs.tsv \
~/ss/allele_Fs_codingRegions_SIFToutput/sorted_GA1_minor_alleles_MISSENSE_DELETERIOUS_HIGH_CONFIDENCE_X.txt | \
awk 'BEGIN { OFS="\t" } { if ($4 == $11) $12 = $5; else $12 = 1 - $5; print }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_DELETERIOUS_X_minor_alleleFs.tsv


# # MISSENSE TOL
join -1 9 -2 6 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.6 2.5 -e "NA" -a 1 -t $'\t' \
${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_TOLERATED_autosomes_alleleFs.tsv \
~/ss/allele_Fs_codingRegions_SIFToutput/sorted_GA1_minor_alleles_MISSENSE_TOLERATED_autosomes.txt | \
awk 'BEGIN { OFS="\t" } { if ($4 == $11) $12 = $5; else $12 = 1 - $5; print }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_TOLERATED_autosomes_minor_alleleFs.tsv

join -1 9 -2 6 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.6 2.5 -e "NA" -a 1 -t $'\t' \
${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_TOLERATED_X_alleleFs.tsv \
~/ss/allele_Fs_codingRegions_SIFToutput/sorted_GA1_minor_alleles_MISSENSE_TOLERATED_X.txt | \
awk 'BEGIN { OFS="\t" } { if ($4 == $11) $12 = $5; else $12 = 1 - $5; print }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_MISSENSE_TOLERATED_X_minor_alleleFs.tsv


# # NONSENSE
join -1 9 -2 6 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.6 2.5 -e "NA" -a 1 -t $'\t' \
${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_NONSENSE_autosomes_alleleFs.tsv \
~/ss/allele_Fs_codingRegions_SIFToutput/sorted_GA1_minor_alleles_NONSENSE_autosomes.txt | \
awk 'BEGIN { OFS="\t" } { if ($4 == $11) $12 = $5; else $12 = 1 - $5; print }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_NONSENSE_autosomes_minor_alleleFs.tsv

join -1 9 -2 6 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.6 2.5 -e "NA" -a 1 -t $'\t' \
${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_NONSENSE_X_alleleFs.tsv \
~/ss/allele_Fs_codingRegions_SIFToutput/sorted_GA1_minor_alleles_NONSENSE_X.txt |\
awk 'BEGIN { OFS="\t" } { if ($4 == $11) $12 = $5; else $12 = 1 - $5; print }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_NONSENSE_X_minor_alleleFs.tsv


# # SYNONYMOUS
join -1 9 -2 6 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.6 2.5 -e "NA" -a 1 -t $'\t' \
${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_SYNONYMOUS_autosomes_alleleFs.tsv \
~/ss/allele_Fs_codingRegions_SIFToutput/sorted_GA1_minor_alleles_SYNONYMOUS_autosomes.txt | \
awk 'BEGIN { OFS="\t" } { if ($4 == $11) $12 = $5; else $12 = 1 - $5; print }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_SYNONYMOUS_autosomes_minor_alleleFs.tsv

join -1 9 -2 6 -o 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.6 2.5 -e "NA" -a 1 -t $'\t' \
${OUTPUT_FOLDER}/sorted_${POPS[${SLURM_ARRAY_TASK_ID}]}_SYNONYMOUS_X_alleleFs.tsv \
~/ss/allele_Fs_codingRegions_SIFToutput/sorted_GA1_minor_alleles_SYNONYMOUS_X.txt | \
awk 'BEGIN { OFS="\t" } { if ($4 == $11) $12 = $5; else $12 = 1 - $5; print }' > ${OUTPUT_FOLDER}/${POPS[${SLURM_ARRAY_TASK_ID}]}_SYNONYMOUS_X_minor_alleleFs.tsv




echo "End of script"

