#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH -t 48:00:00                               # Set time limit to 36 hours
#SBATCH --job-name=ss.genomewide_heterozygosity_PLINK        #Job name
#SBATCH --mem 20G
#SBATCH -o ss.genomewide_heterozygosity_PLINK.out                    #Standard output log
#SBATCH -e ss.genomewide_heterozygosity_PLINK.err                     #Standard error log


# load modules
module load vcftools/0.1.16
module load plink/1.90

MASTER=~/ss
VCF=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic.vcf.gz
VCF_onlyLGAUTOSOMES=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_autosomes.vcf.gz

mkdir -p ${MASTER}/data/data_for_plink

# # filter VCF to autosomes
# # commented out because I've already made this file for another analysis
# vcftools \
# --gzvcf ${VCF} \
# --chr LG10 \
# --chr LG2 \
# --chr LG3 \
# --chr LG4 \
# --chr LG5 \
# --chr LG6 \
# --chr LG7 \
# --chr LG8 \
# --chr LG9 \
# --recode \
# --recode-INFO-all \
# --out ${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyLGautosomes

# gzip ${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyLGautosomes.recode.vcf
# mv ${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyLGautosomes.recode.vcf.gz \
# ${VCF_onlyLGAUTOSOMES}

#convert VCF to .bed
plink \
--allow-extra-chr \
--set-missing-var-ids @:# \
--vcf ${VCF_onlyLGAUTOSOMES} \
--make-bed \
--out ${MASTER}/data/data_for_plink/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyLGautosomes # plink will add a suffix to this



# compuate the heterozygosity
# you provide plink the base_name of the 3 files produced above, and it looks for each of the .bed, .bim, .fam files
plink \
--allow-extra-chr \
--set-missing-var-ids @:# \
--bfile ${MASTER}/data/data_for_plink/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyLGautosomes \
--het \
--out ${MASTER}/plink_het/ss_genomewide_heterozygosity




# output:

# FID IID O(HET) O(HOM) N_NM
# 1   1   0.122   0.878   1000
# 1   2   0.134   0.866   1000
# ...
# family ID
# indiv ID
# Observed Het (proportion)
# Observed Hom (proportion)
# Number non-missing genptypes

