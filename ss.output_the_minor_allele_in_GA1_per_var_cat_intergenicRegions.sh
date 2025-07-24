#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 10G
##SBATCH --cpus-per-task 4
#SBATCH -t 150:00:00                               # Set time limit 
#SBATCH --job-name=ss.output_the_minor_allele_in_GA1_per_var_cat_intergenicRegions       #Job name
#SBATCH -o ss.output_the_minor_allele_in_GA1_per_var_cat_intergenicRegions.out           #Standard output log
#SBATCH -e ss.output_the_minor_allele_in_GA1_per_var_cat_intergenicRegions.err           #Standard error log


# # set up the environment
module add bcftools
module add htslib/1.15.1


MASTER=~/ss

INTERGENIC_SNP_VCF=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_GA1_mono_poly__intergenicRegions.vcf.gz

OUTPUT_FOLDER=${MASTER}/allele_Fs_intergenicRegions

PATH_TO_POPFILES=${MASTER}/population_files/182sample_popfiles
GA1_POP_FILE=${PATH_TO_POPFILES}/ss.samples_GA1.txt


### for each variant in a VCF, find the minor allele in the GA1 popualtion at that locus

bcftools view -S ${GA1_POP_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${INTERGENIC_SNP_VCF} | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
awk 'BEGIN { OFS="\t" } { if ($5 <= 0.5) print $1, $2, $3, $4, $4, $1"_"$2; else print $1, $2, $3, $4, $3, $1"_"$2 }' > ${OUTPUT_FOLDER}/GA1_minor_alleles_INTERGENIC_autosomes.txt
sort -k6,6 ${OUTPUT_FOLDER}/GA1_minor_alleles_INTERGENIC_autosomes.txt > ${OUTPUT_FOLDER}/sorted_GA1_minor_alleles_INTERGENIC_autosomes.txt

bcftools view -S ${GA1_POP_FILE} -r LGX ${INTERGENIC_SNP_VCF} | \
bcftools +fill-tags -- -t AF | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n' | \
awk 'BEGIN { OFS="\t" } { if ($5 <= 0.5) print $1, $2, $3, $4, $4, $1"_"$2; else print $1, $2, $3, $4, $3, $1"_"$2 }' > ${OUTPUT_FOLDER}/GA1_minor_alleles_INTERGENIC_X.txt
sort -k6,6 ${OUTPUT_FOLDER}/GA1_minor_alleles_INTERGENIC_X.txt > ${OUTPUT_FOLDER}/sorted_GA1_minor_alleles_INTERGENIC_X.txt






# remove unsorted files
rm \
${OUTPUT_FOLDER}/GA1_minor_alleles_INTERGENIC_autosomes.txt \
${OUTPUT_FOLDER}/GA1_minor_alleles_INTERGENIC_X.txt