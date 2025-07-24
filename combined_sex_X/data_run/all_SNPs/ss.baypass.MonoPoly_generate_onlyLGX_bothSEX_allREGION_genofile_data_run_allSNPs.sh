#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH -t 48:00:00                               # Set time limit to 36 hours
#SBATCH --job-name=ss.baypass.MonoPoly_bothSEX_generate_onlyLGX_allREGION_genofile_data_run_allSNPs     #Job name
#SBATCH --mem 30G
#SBATCH -o ss.baypass.MonoPoly_bothSEX_generate_onlyLGX_allREGION_genotfile_data_run_allSNPs.out                    #Standard output log
#SBATCH -e ss.baypass.MonoPoly_bothSEX_generate_onlyLGX_allREGION_genotfile_data_run_allSNPs.err                     #Standard error log



#set up environment
module add vcftools/0.1.16
module add plink/1.90
module add bcftools/1.15.1

MASTER=~/ss
VCF=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic.vcf.gz
MONO_POLY_SAMPLES=${MASTER}/population_files/182sample_popfiles/ss.samples_mono_poly.txt

mkdir -p ${MASTER}/data/vcfs_for_data_run_MonoPoly_bothSEX_onlyLGX_allREGION_allSNPs
WORKING_FOLDER=${MASTER}/data/vcfs_for_data_run_MonoPoly_bothSEX_onlyLGX_allREGION_allSNPs





##################################################### COUNT SNPS IN RAW VCF TO COMPARE TO POST-FILTER
N_SNPS=$(zgrep -vc "^#" "${VCF}")
echo "RAW VCF N_snps=${N_SNPS}"




##################################################### FILTER TO Mono / Poly SS TREATMENTS
##
## INPUT:  ANNOTAED RAW VCF
## OUTPUT: onlyMonoPoly
vcftools \
--gzvcf ${VCF} \
--keep ${MONO_POLY_SAMPLES} \
--non-ref-ac-any 1 \
--recode \
--recode-INFO-all \
--out ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX


# gzip output
gzip ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX.recode.vcf
#rename file
mv ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX.recode.vcf.gz \
${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX.vcf.gz


N_SNPS=$(zgrep -vc "^#" "${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX.vcf.gz")
echo "vcftools raw VCF filter to only Mono & Poly MALE samples & remove invariant sites - DONE... N_snps=${N_SNPS}"

VCF_allREGION_onlyMonoPoly_bothSEX=${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX.vcf.gz





##################################################### FILTER TO CHROMOSOME-LEVEL SCAFFOLDS
##
## INPUT:  VCF onlyMonoPoly_bothSEX
## OUTPUT: VCF onlyMonoPoly_bothSEX onlyLGX
vcftools \
--gzvcf ${VCF_allREGION_onlyMonoPoly_bothSEX} \
--chr LGX \
--recode \
--recode-INFO-all \
--out ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX

# gzip output
gzip ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX.recode.vcf
mv  ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX.recode.vcf.gz \
${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX.vcf.gz

N_SNPS=$(zgrep -vc "^#" "${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX.vcf.gz")
echo "vcftools allREGION onlyMonoPoly_bothSEX filter to onlyLGX - DONE... N_snps=${N_SNPS}"

VCF_allREGION_onlyMonoPoly_bothSEX_onlyLGX=${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX.vcf.gz





##################################################### ADD INFORMATION TO VCF FILE
##
## INPUT:  RAW VCF
## OUTPUT: ANNOTATED RAW VCF
bcftools \
+fill-tags \
${VCF_allREGION_onlyMonoPoly_bothSEX_onlyLGX} \
-Oz -o ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX_TAGS.vcf.gz \
-- -t MAF

N_SNPS=$(zgrep -vc "^#" "${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX_TAGS.vcf.gz")
echo "bcftools add tags to raw VCF - DONE... N_snps=${N_SNPS}"

VCF_allREGION_onlyMonoPoly_bothSEX_onlyLGX_TAGS=${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX_TAGS.vcf.gz





##################################################### WHOLE DATASET MAF FILTER
##
## INPUT:  VCF onlyLGX
## OUTPUT: VCF onlyH_L onlyLGX MAF>0.03
bcftools \
view -e 'MAF<0.04' \
-Oz --output ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX_MAF003.vcf.gz \
${VCF_allREGION_onlyMonoPoly_bothSEX_onlyLGX_TAGS}

bcftools index ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX_MAF003.vcf.gz

N_SNPS=$(zgrep -vc "^#" "${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX_MAF003.vcf.gz")
echo "bcftools onlyMonoPoly_bothSEX onlyLGX filter to MAF>0.03 - DONE... N_snps=${N_SNPS}"

VCF_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF=${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX_MAF003.vcf.gz





# ##################################################### WHOLE DATASET LD PRUNE
# ##
# ## INPUT: VCF onlyH_L onlyLGX MAF0.03
# ## OUTPUT: VCF onlyH_L onlyLGX MAF0.03 LDPRUNE_50KB_0.3
# bcftools \
# +prune -m 0.3 -w 50kb \
# -Oz -o ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX_MAF003_LD3050.vcf.gz \
# ${VCF_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF}

# N_SNPS=$(zgrep -vc "^#" "${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX_MAF003_LD3050.vcf.gz")
# echo "bcftools onlyMonoPoly_bothSEX onlyLGX MAF0.03 LDprune 50kb 0.3 - DONE... N_snps=${N_SNPS}"

# VCF_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_LD0350=${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_onlyMonoPoly_bothSEX_onlyLGX_MAF003_LD3050.vcf.gz





# extract the positions from this file
mkdir -p ${MASTER}/baypass/positions_of_snps_used_baypass_data_runs/
zgrep -v "^#" ${VCF_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF} |
awk -F'\t' '{print $1"\t"$2}' > \
${MASTER}/baypass/positions_of_snps_used_baypass_data_runs/allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune.txt
echo "Extract the snp positions in this final filtered vcf and save to file - DONE."



MONO_POLY_POPS=(monoA monoB monoC polyA polyB polyC)
#
for POP in ${MONO_POLY_POPS[@]}
do
    
    # take VCF after all previous filters
    # filter samples to a single focal population
    # extract the genotype counts
    vcftools \
    --gzvcf ${VCF_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF} \
    --keep ${MASTER}/population_files/182sample_popfiles/ss.samples_${POP}.txt \
    --counts2 \
    --out ${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_${POP}_genotypes_tmp

    # get the relevant columns and change the delim to a space
    tail -n+2 ${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_${POP}_genotypes_tmp.frq.count | \
    cut -f 5,6 --output-delimiter=" " > \
    ${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_${POP}_geno

    # remove the file with the wrong delimiter
    rm ${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_${POP}_genotypes_tmp.frq.count

done


# PASTE GENO FILES FOR EACH POP BACK TOGETHER
paste \
${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_monoA_geno \
${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_monoB_geno \
${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_monoC_geno \
${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_polyA_geno \
${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_polyB_geno \
${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_polyC_geno -d" " > \
${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_genotypes.txt

rm -f ${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_?????_geno

N_GENOTYPES=$(cat "${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_genotypes.txt" | wc -l)
N_COLS=$(head -n 1 "${MASTER}/data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_genotypes.txt" | wc -w)
echo "The final genotype file data/genofiles_for_data_run/6pop_allREGION_onlyMonoPoly_bothSEX_onlyLGX_MAF003_noLDprune_genotypes.txt contains ${N_GENOTYPES} rows and ${N_COLS} columns"


echo "end of script"