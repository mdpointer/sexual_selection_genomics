#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 10G
##SBATCH --cpus-per-task 4
#SBATCH -t 150:00:00                               # Set time limit
#SBATCH --job-name=ss.output_the_minor_allele_in_GA1_per_var_cat_codingRegions       #Job name
#SBATCH -o ss.output_the_minor_allele_in_GA1_per_var_cat_codingRegions.out           #Standard output log
#SBATCH -e ss.output_the_minor_allele_in_GA1_per_var_cat_codingRegions.err           #Standard error log


# # set up the environment
module add bcftools
module add htslib/1.15.1


MASTER=~/ss

MISSENSE_DELETERIOUS_VCF_HIGH_CONFIDENCE=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_SIFT_MISSENSE_DELETERIOUS_high_confidence.vcf.gz
MISSENSE_TOLERATED_VCF=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_SIFT_MISSENSE_TOLERATED.vcf.gz
NONSENSE_SNP_VCF=${MASTER}/snpEff/codingRegions_snpEff_NONSENSE.vcf.gz
SYNONYMOUS_SNP_VCF=${MASTER}/snpEff/codingRegions_snpEff_SYNONYMOUS.vcf.gz

OUTPUT_FOLDER=${MASTER}/alelle_Fs_codingRegions_SIFToutput

PATH_TO_POPFILES=${MASTER}/population_files/182sample_popfiles
GA1_POP_FILE=${PATH_TO_POPFILES}/ss.samples_GA1.txt


### for each variant category, find the minor allele in the GA1 popualtion at each locus

# MISSENSE DEL
bcftools view -S ${GA1_POP_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${MISSENSE_DELETERIOUS_VCF_HIGH_CONFIDENCE} | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
awk 'BEGIN { OFS="\t" } { if ($5 <= 0.5) print $1, $2, $3, $4, $4, $1"_"$2; else print $1, $2, $3, $4, $3, $1"_"$2 }' > ${OUTPUT_FOLDER}/GA1_minor_alleles_MISSENSE_DELETERIOUS_HIGH_CONFIDENCE_autosomes.txt
sort -k6,6 ${OUTPUT_FOLDER}/GA1_minor_alleles_MISSENSE_DELETERIOUS_HIGH_CONFIDENCE_autosomes.txt > ${OUTPUT_FOLDER}/sorted_GA1_minor_alleles_MISSENSE_DELETERIOUS_HIGH_CONFIDENCE_autosomes.txt

bcftools view -S ${GA1_POP_FILE} -r LGX ${MISSENSE_DELETERIOUS_VCF_HIGH_CONFIDENCE} | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
awk 'BEGIN { OFS="\t" } { if ($5 <= 0.5) print $1, $2, $3, $4, $4, $1"_"$2; else print $1, $2, $3, $4, $3, $1"_"$2 }' > ${OUTPUT_FOLDER}/GA1_minor_alleles_MISSENSE_DELETERIOUS_HIGH_CONFIDENCE_X.txt
sort -k6,6 ${OUTPUT_FOLDER}/GA1_minor_alleles_MISSENSE_DELETERIOUS_HIGH_CONFIDENCE_X.txt > ${OUTPUT_FOLDER}/sorted_GA1_minor_alleles_MISSENSE_DELETERIOUS_HIGH_CONFIDENCE_X.txt


# MISSENSE TOL
bcftools view -S ${GA1_POP_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${MISSENSE_TOLERATED_VCF} | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
awk 'BEGIN { OFS="\t" } { if ($5 <= 0.5) print $1, $2, $3, $4, $4, $1"_"$2; else print $1, $2, $3, $4, $3, $1"_"$2 }' > ${OUTPUT_FOLDER}/GA1_minor_alleles_MISSENSE_TOLERATED_autosomes.txt
sort -k6,6 ${OUTPUT_FOLDER}/GA1_minor_alleles_MISSENSE_TOLERATED_autosomes.txt > ${OUTPUT_FOLDER}/sorted_GA1_minor_alleles_MISSENSE_TOLERATED_autosomes.txt

bcftools view -S ${GA1_POP_FILE} -r LGX ${MISSENSE_TOLERATED_VCF} | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
awk 'BEGIN { OFS="\t" } { if ($5 <= 0.5) print $1, $2, $3, $4, $4, $1"_"$2; else print $1, $2, $3, $4, $3, $1"_"$2 }' > ${OUTPUT_FOLDER}/GA1_minor_alleles_MISSENSE_TOLERATED_X.txt
sort -k6,6 -k2,2n ${OUTPUT_FOLDER}/GA1_minor_alleles_MISSENSE_TOLERATED_X.txt > ${OUTPUT_FOLDER}/sorted_GA1_minor_alleles_MISSENSE_TOLERATED_X.txt


# NONSENSE
bcftools view -S ${GA1_POP_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${NONSENSE_SNP_VCF} | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
awk 'BEGIN { OFS="\t" } { if ($5 <= 0.5) print $1, $2, $3, $4, $4, $1"_"$2; else print $1, $2, $3, $4, $3, $1"_"$2 }' > ${OUTPUT_FOLDER}/GA1_minor_alleles_NONSENSE_autosomes.txt
sort -k6,6 ${OUTPUT_FOLDER}/GA1_minor_alleles_NONSENSE_autosomes.txt > ${OUTPUT_FOLDER}/sorted_GA1_minor_alleles_NONSENSE_autosomes.txt

bcftools view -S ${GA1_POP_FILE} -r LGX ${NONSENSE_SNP_VCF} | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
awk 'BEGIN { OFS="\t" } { if ($5 <= 0.5) print $1, $2, $3, $4, $4, $1"_"$2; else print $1, $2, $3, $4, $3, $1"_"$2 }' > ${OUTPUT_FOLDER}/GA1_minor_alleles_NONSENSE_X.txt
sort -k6,6 ${OUTPUT_FOLDER}/GA1_minor_alleles_NONSENSE_X.txt > ${OUTPUT_FOLDER}/sorted_GA1_minor_alleles_NONSENSE_X.txt


# SYNONYMOUS
bcftools view -S ${GA1_POP_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${SYNONYMOUS_SNP_VCF} | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
awk 'BEGIN { OFS="\t" } { if ($5 <= 0.5) print $1, $2, $3, $4, $4, $1"_"$2; else print $1, $2, $3, $4, $3, $1"_"$2 }' > ${OUTPUT_FOLDER}/GA1_minor_alleles_SYNONYMOUS_autosomes.txt
sort -k6,6 ${OUTPUT_FOLDER}/GA1_minor_alleles_SYNONYMOUS_autosomes.txt > ${OUTPUT_FOLDER}/sorted_GA1_minor_alleles_SYNONYMOUS_autosomes.txt

bcftools view -S ${GA1_POP_FILE} -r LGX ${SYNONYMOUS_SNP_VCF} | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
awk 'BEGIN { OFS="\t" } { if ($5 <= 0.5) print $1, $2, $3, $4, $4, $1"_"$2; else print $1, $2, $3, $4, $3, $1"_"$2 }' > ${OUTPUT_FOLDER}/GA1_minor_alleles_SYNONYMOUS_X.txt
sort -k6,6 ${OUTPUT_FOLDER}/GA1_minor_alleles_SYNONYMOUS_X.txt > ${OUTPUT_FOLDER}/sorted_GA1_minor_alleles_SYNONYMOUS_X.txt



# remove unsorted files
rm \
${OUTPUT_FOLDER}/GA1_minor_alleles_MISSENSE_DELETERIOUS_HIGH_CONFIDENCE_autosomes.txt \
${OUTPUT_FOLDER}/GA1_minor_alleles_MISSENSE_DELETERIOUS_HIGH_CONFIDENCE_X.txt \
${OUTPUT_FOLDER}/GA1_minor_alleles_MISSENSE_TOLERATED_autosomes.txt \
${OUTPUT_FOLDER}/GA1_minor_alleles_MISSENSE_TOLERATED_X.txt \
${OUTPUT_FOLDER}/GA1_minor_alleles_NONSENSE_autosomes.txt \
${OUTPUT_FOLDER}/GA1_minor_alleles_NONSENSE_X.txt \
${OUTPUT_FOLDER}/GA1_minor_alleles_SYNONYMOUS_autosomes.txt \
${OUTPUT_FOLDER}/GA1_minor_alleles_SYNONYMOUS_X.txt