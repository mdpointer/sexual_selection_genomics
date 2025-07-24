#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH -t 48:00:00                               # Set time limit to 36 hours
#SBATCH --job-name=ss.baypass.MonoPoly_generate_onlyLGautosomes_nonEXONIC10kb_genotype_matrix        #Job name
#SBATCH --mem 30G
#SBATCH -o ss.baypass.MonoPoly_generate_onlyLGautosomes_nonEXONIC10kb_genotype_matrix.out                    #Standard output log
#SBATCH -e ss.baypass.MonoPoly_generate_onlyLGautosomes_nonEXONIC10kb_genotype_matrix.err                     #Standard error log



#set up environment
module add vcftools/0.1.16
module add plink/1.90
module add bcftools/1.15.1

MASTER=~/ss
VCF=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC.vcf.gz
mono_poly_SAMPLES=${MASTER}/population_files/182sample_popfiles/ss.samples_mono_poly.txt

mkdir -p ${MASTER}/data/vcfs_for_omega_run_MonoPoly_onlyLGautosomes_nonEXONIC10kb
WORKING_FOLDER=${MASTER}/data/vcfs_for_omega_run_MonoPoly_onlyLGautosomes_nonEXONIC10kb


##################################################### COUNT SNPS IN RAW VCF TO COMPARE TO POST-FILTER
N_SNPS=$(zgrep -vc "^#" "${VCF}")
echo "RAW VCF N_snps=${N_SNPS}"




##################################################### FILTER TO MONO / POLY DISPERSAL TREATMENTS
##
## INPUT:  ANNOTAED RAW VCF
## OUTPUT: onlyMonoPoly
vcftools \
--gzvcf ${VCF} \
--keep ${mono_poly_SAMPLES} \
--non-ref-ac-any 1 \
--recode \
--recode-INFO-all \
--out ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly


# gzip output
gzip ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly.recode.vcf
#rename file
mv ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly.recode.vcf.gz \
${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly.vcf.gz


N_SNPS=$(zgrep -vc "^#" "${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly.vcf.gz")
echo "vcftools raw VCF filter to only Mono & Poly samples & remove invariant sites - DONE... N_snps=${N_SNPS}"

VCF_nonEXONIC10kb_onlyMP=${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly.vcf.gz





##################################################### FILTER TO CHROMOSOME-LEVEL SCAFFOLDS
##
## INPUT:  VCF onlyMonoPoly
## OUTPUT: VCF onlyMonoPoly onlyLGautosomes
vcftools \
--gzvcf ${VCF_nonEXONIC10kb_onlyMP} \
--chr LG10 \
--chr LG2 \
--chr LG3 \
--chr LG4 \
--chr LG5 \
--chr LG6 \
--chr LG7 \
--chr LG8 \
--chr LG9 \
--recode \
--recode-INFO-all \
--out ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes

# gzip output
gzip ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes.recode.vcf
mv  ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes.recode.vcf.gz \
${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes.vcf.gz

N_SNPS=$(zgrep -vc "^#" "${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes.vcf.gz")
echo "vcftools nonEXONIC10kb onlyMonoPoly filter to onlyLGautosomes - DONE... N_snps=${N_SNPS}"

VCF_nonEXONIC10kb_onlyMP_onlyLGautosomes=${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes.vcf.gz





##################################################### ADD INFORMATION TO VCF FILE
##
## INPUT:  RAW VCF
## OUTPUT: ANNOTATED RAW VCF
bcftools \
+fill-tags \
${VCF_nonEXONIC10kb_onlyMP_onlyLGautosomes} \
-Oz -o ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_tags.vcf.gz \
-- -t MAF

N_SNPS=$(zgrep -vc "^#" "${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_tags.vcf.gz")
echo "bcftools add tags to raw VCF - DONE... N_snps=${N_SNPS}"

VCF_nonEXONIC10kb_onlyMP_onlyLGautosomes_TAGS=${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_tags.vcf.gz





##################################################### WHOLE DATASET MAF FILTER
##
## INPUT:  VCF onlyLGautosomes
## OUTPUT: VCF onlyH_L onlyLGautosomes MAF>0.03
bcftools \
view -e 'MAF<0.04' \
-Oz --output ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003.vcf.gz \
${VCF_nonEXONIC10kb_onlyMP_onlyLGautosomes_TAGS}

bcftools index ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003.vcf.gz

N_SNPS=$(zgrep -vc "^#" "${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003.vcf.gz")
echo "bcftools onlyMonoPoly onlyLGautosomes filter to MAF>0.03 - DONE... N_snps=${N_SNPS}"

VCF_nonEXONIC10kb_onlyMP_onlyLGautosomes_MAF=${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003.vcf.gz





##################################################### REMOVE SITES IN <6 INDS IN EVERY POPULATION
##
## INPUT: VCF onlyH_L onlyLGautosomes MAF0.03 LDPRUNE_50KB_0.3
## OUTPUT: VCF onlyH_L onlyLGautosomes MAF0.03 LDPRUNE_50KB_0.3 SNPs_in_>=4INDSperPOP

MONO_POLY_POPS=(monoA monoB monoC polyA polyB polyC)
mkdir -p ${WORKING_FOLDER}/indiv_pop_vcfs
mkdir -p ${MASTER}/data/genofiles_for_omega_run

## SPLIT INTO 1 VCF PER POP (SOME SORT OF LOOP, BCF/VCFTOOLS)
for POP in ${MONO_POLY_POPS[@]}
do
    # output population-specific VCF
    bcftools view \
    -S ${MASTER}/population_files/182sample_popfiles/ss.samples_${POP}.txt \
    -Oz -o ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_${POP}.vcf.gz \
    ${VCF_nonEXONIC10kb_onlyMP_onlyLGautosomes_MAF}

    # remove sites missing in more than 2 individuals in the focal pop
    vcftools \
    --gzvcf ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_${POP}.vcf.gz \
    --max-missing-count 5 \
    --recode \
    --recode-INFO-all \
    --out ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_${POP}_MAXMISS5


    # gzip
    gzip ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_${POP}_MAXMISS5.recode.vcf
    # rename
    mv ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_${POP}_MAXMISS5.recode.vcf.gz \
    ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_${POP}_MAXMISS5.vcf.gz

    # extract positions of SNPs remaining in the focal populatio VCF
    zgrep -v "^#" ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_${POP}_MAXMISS5.vcf.gz | \
    awk '{print $1 "\t" $2}' | \
    sort > \
    ${WORKING_FOLDER}/indiv_pop_vcfs/SNPs.nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_${POP}_MAXMISS5_SNPpos.txt

    # echo the number of snps remaining in the focal population after the MAX MISSING filter
    N_SNPS=$(grep -vc "^#" "${MASTER}/data/vcfs_for_omega_run_MonoPoly_onlyLGautosomes_nonEXONIC10kb/indiv_pop_vcfs/SNPs.nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_${POP}_MAXMISS5_SNPpos.txt")
    echo "MAX MISSING FILTER IN ${POP} - DONE... N_snps=${N_SNPS}"

    # remove the MAXMIS VCF for focal pop, once the snp positions have been recorded
    rm -f ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_${POP}_MAXMISS5.vcf.gz

done




# TAKE THE OUTPUT FILES FROM ABOVE AND FIND POSITIONS COMMON TO ALL POPULATIONS (to use comm on >2 files requires temp files)
mkdir -p ${MASTER}/tmp

comm -12 ${WORKING_FOLDER}/indiv_pop_vcfs/SNPs.nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_${MONO_POLY_POPS[0]}_MAXMISS5_SNPpos.txt ${WORKING_FOLDER}/indiv_pop_vcfs/SNPs.nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_${MONO_POLY_POPS[1]}_MAXMISS5_SNPpos.txt > ${MASTER}/tmp/comm_temp01.txt
comm -12 ${MASTER}/tmp/comm_temp01.txt ${WORKING_FOLDER}/indiv_pop_vcfs/SNPs.nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_${MONO_POLY_POPS[2]}_MAXMISS5_SNPpos.txt > ${MASTER}/tmp/comm_temp02.txt
comm -12 ${MASTER}/tmp/comm_temp02.txt ${WORKING_FOLDER}/indiv_pop_vcfs/SNPs.nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_${MONO_POLY_POPS[3]}_MAXMISS5_SNPpos.txt > ${MASTER}/tmp/comm_temp03.txt
comm -12 ${MASTER}/tmp/comm_temp03.txt ${WORKING_FOLDER}/indiv_pop_vcfs/SNPs.nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_${MONO_POLY_POPS[4]}_MAXMISS5_SNPpos.txt > ${MASTER}/tmp/comm_temp04.txt
comm -12 ${MASTER}/tmp/comm_temp04.txt ${WORKING_FOLDER}/indiv_pop_vcfs/SNPs.nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_${MONO_POLY_POPS[5]}_MAXMISS5_SNPpos.txt > ${WORKING_FOLDER}/snps_common_to_all_6pops_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5.txt



nonEXONIC10kb_onlyLGautosomes_SNPS_PRESENT_IN_ALL_POPS=${WORKING_FOLDER}/snps_common_to_all_6pops_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5.txt

# ECHO THE NUMBER OF SNPS COMMON TO ALL POPULATIONS AFTER THE MAX MISSING FITLER HAS BEEN APPLIED TO EACH POPULATION INDIVIDUALLY
N_SNPS=$(grep -vc "^#" "${WORKING_FOLDER}/snps_common_to_all_6pops_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5.txt")
echo "MAX MISSING IN EACH POP & FIND onlyLGautosome SNPs COMMON TO ALL MonoPoly POPS - DONE... N_snps=${N_SNPS}"

# remove all temporary files
rm -f ${MASTER}/tmp/comm_temp??.txt



# FILTER THE PRE-MAXMISSING VCF TO SNPS COMMON TO ALL POPS AFTER THE MAXMISS APPLIED
bcftools \
view \
-R ${nonEXONIC10kb_onlyLGautosomes_SNPS_PRESENT_IN_ALL_POPS} \
-Oz -o ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5.vcf.gz \
${VCF_nonEXONIC10kb_onlyMP_onlyLGautosomes_MAF}

N_SNPS=$(zgrep -vc "^#" "${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5.vcf.gz")
echo "bcftools onlyMonoPoly onlyLGautosomes MAF0.03 MAXMISS5 - DONE... N_snps=${N_SNPS}"

VCF_nonEXONIC10kb_onlyMP_onlyLGautosomes_MAF003_MAXMISS5=${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5.vcf.gz




##################################################### WHOLE DATASET LD PRUNE
##
## INPUT: VCF onlyH_L onlyLGautosomes MAF0.03
## OUTPUT: VCF onlyH_L onlyLGautosomes MAF0.03 LDPRUNE_50KB_0.3
bcftools \
+prune -m 0.3 -w 50kb \
-Oz -o ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD3050.vcf.gz \
${VCF_nonEXONIC10kb_onlyMP_onlyLGautosomes_MAF003_MAXMISS5}

N_SNPS=$(zgrep -vc "^#" "${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD3050.vcf.gz")
echo "bcftools onlyMonoPoly onlyLGautosomes MAF0.03 MAXISS5 LDprune 50kb 0.3 - DONE... N_snps=${N_SNPS}"

VCF_nonEXONIC10kb_onlyMP_onlyLGautosomes_MAF003_MAXMISS5_LD0350=${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD3050.vcf.gz





# extract the positions from this file
mkdir -p ${MASTER}/baypass/positions_of_snps_used_generate_omega_matrices/
zgrep -v "^#" ${WORKING_FOLDER}/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_NONEXONIC_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD3050.vcf.gz |
awk -F'\t' '{print $1"\t"$2}' > \
${MASTER}/baypass/positions_of_snps_used_generate_omega_matrices/nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350.txt
echo "Extract the snp positions in this final filtered vcf and save to file - DONE."




#
for POP in ${MONO_POLY_POPS[@]}
do
    
    # take VCF after all previous filters
    # filter samples to a single focal population
    # extract the genotype counts
    vcftools \
    --gzvcf ${VCF_nonEXONIC10kb_onlyMP_onlyLGautosomes_MAF003_MAXMISS5_LD0350} \
    --keep ${MASTER}/population_files/182sample_popfiles/ss.samples_${POP}.txt \
    --counts2 \
    --out ${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_${POP}_genotypes_tmp

    # get the relevant columns and change the delim to a space
    tail -n+2 ${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_${POP}_genotypes_tmp.frq.count | \
    cut -f 5,6 --output-delimiter=" " > \
    ${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_${POP}_geno

    # remove the file with the wrong delimiter
    rm ${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_${POP}_genotypes_tmp.frq.count

done


# PASTE GENO FILES FOR EACH POP BACK TOGETHER
paste \
${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_monoA_geno \
${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_monoB_geno \
${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_monoC_geno \
${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_polyA_geno \
${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_polyB_geno \
${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_polyC_geno -d" " > \
${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt

rm -f ${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_?????_geno

N_GENOTYPES=$(cat "${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt" | wc -l)
N_COLS=$(head -n 1 "${MASTER}/data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt" | wc -w)
echo "The final genotype file data/genofiles_for_omega_run/6pop_nonEXONIC10kb_onlyMonoPoly_onlyLGautosomes_MAF003_MAXMISS5_LD0350_genotypes.txt contains ${N_GENOTYPES} rows and ${N_COLS} columns"


echo "end of script"