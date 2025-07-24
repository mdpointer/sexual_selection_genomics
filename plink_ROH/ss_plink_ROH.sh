#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH -t 48:00:00                               # Set time limit to 36 hours
#SBATCH --job-name=ss.ROH_PLINK        #Job name
#SBATCH --mem 20G
#SBATCH -o ss.ROH_PLINK.out                    #Standard output log
#SBATCH -e ss.ROH_PLINK.err                     #Standard error log


# load modules
module load vcftools/0.1.16
module load plink/1.90


MASTER=~/ss
VCF=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic.vcf.gz
VCF_onlyLGAUTOSOMES=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_autosomes.vcf.gz
VCF_onlyLGAUTOSOMES_unZIPPED=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_autosomes.vcf
BED_onlyLGAUTOSOMES=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_autosomes

OUT_DIR=${MASTER}/plink_ROH

mkdir -p ${MASTER}/data/data_for_plink

mkdir -p ${OUT_DIR}

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


# # unzip the input VCF
# gunzip -c ${VCF_onlyLGAUTOSOMES} > ${VCF_onlyLGAUTOSOMES_unZIPPED}

# #convert VCF to .bed
# plink \
# --allow-extra-chr \
# --set-missing-var-ids @:# \
# --vcf ${VCF_onlyLGAUTOSOMES_unZIPPED} \
# --make-bed \
# --out ${BED_onlyLGAUTOSOMES} # plink will add a suffix to this


# Run plink ROH
plink \
  --bfile ${BED_onlyLGAUTOSOMES} \
  --allow-extra-chr \
  --double-id \
  --chr-set 50 \
  --homozyg \
  --homozyg-snp 30 \
  --homozyg-kb 200 \
  --homozyg-density 100 \
  --homozyg-gap 100 \
  --homozyg-window-snp 30 \
  --homozyg-window-het 2 \
  --homozyg-window-missing 3 \
  --homozyg-window-threshold 0.05 \
  --homozyg-group \
  --out ${OUT_DIR}/plink_ROH_run3



# # remove the uncompressed VCF created above
# rm ${VCF_onlyLGAUTOSOMES_unZIPPED}