#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH -t 48:00:00                               # Set time limit to 36 hours
#SBATCH --job-name=ss.genomewide_diversity_vcftools        #Job name
#SBATCH --mem 20G
#SBATCH -o ss.genomewide_diversity_vcftools.out                    #Standard output log
#SBATCH -e ss.genomewide_diversity_vcftools.err                     #Standard error log



# load modules
module load vcftools/0.1.16

MASTER=~/ss
VCF_onlyLGAUTOSOMES=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_biallelic_autosomes.vcf.gz

mkdir -p ${MASTER}/vcftools_diversity


GA1_popfile=${MASTER}/population_files/182sample_popfiles/ss.samples_GA1.txt
monoA_popfile=${MASTER}/population_files/182sample_popfiles/ss.samples_monoA.txt
monoB_popfile=${MASTER}/population_files/182sample_popfiles/ss.samples_monoB.txt
monoC_popfile=${MASTER}/population_files/182sample_popfiles/ss.samples_monoC.txt
polyA_popfile=${MASTER}/population_files/182sample_popfiles/ss.samples_polyA.txt
polyB_popfile=${MASTER}/population_files/182sample_popfiles/ss.samples_polyB.txt
polyC_popfile=${MASTER}/population_files/182sample_popfiles/ss.samples_polyC.txt

OUT_PATH=${MASTER}/vcftools_diversity
OUT_PREFIX=vcftools_perPop_perSite_pi





# Calculate for GA1
# vcftools --gzvcf ${VCF_onlyLGAUTOSOMES} --keep ${GA1_popfile} --site-pi --out ${OUT_PATH}/${OUT_PREFIX}_GA1
PI_GA1=$(awk 'NR>1 && $3 !~ /-?nan/ && $3+0==$3 {sum+=$3; count++} END {print (count>0 ? sum/count : "NA")}' ${OUT_PATH}/${OUT_PREFIX}_GA1.sites.pi)

# Calculate for monoA
# vcftools --gzvcf ${VCF_onlyLGAUTOSOMES} --keep "$monoA_popfile" --site-pi --out ${OUT_PATH}/${OUT_PREFIX}_monoA
PI_MONOA=$(awk 'NR>1 && $3 !~ /-?nan/ && $3+0==$3 {sum+=$3; count++} END {print (count>0 ? sum/count : "NA")}' ${OUT_PATH}/${OUT_PREFIX}_monoA.sites.pi)

# Calculate for monoB
# vcftools --gzvcf ${VCF_onlyLGAUTOSOMES} --keep "$monoB_popfile" --site-pi --out ${OUT_PATH}/${OUT_PREFIX}_monoB
PI_MONOB=$(awk 'NR>1 && $3 !~ /-?nan/ && $3+0==$3 {sum+=$3; count++} END {print (count>0 ? sum/count : "NA")}' ${OUT_PATH}/${OUT_PREFIX}_monoB.sites.pi)

# Calculate for monoC
# vcftools --gzvcf ${VCF_onlyLGAUTOSOMES} --keep "$monoC_popfile" --site-pi --out ${OUT_PATH}/${OUT_PREFIX}_monoC
PI_MONOC=$(awk 'NR>1 && $3 !~ /-?nan/ && $3+0==$3 {sum+=$3; count++} END {print (count>0 ? sum/count : "NA")}' ${OUT_PATH}/${OUT_PREFIX}_monoC.sites.pi)

# Calculate for polyA
# vcftools --gzvcf ${VCF_onlyLGAUTOSOMES} --keep "$polyA_popfile" --site-pi --out ${OUT_PATH}/${OUT_PREFIX}_polyA
PI_POLYA=$(awk 'NR>1 && $3 !~ /-?nan/ && $3+0==$3 && $3 !~ /-?nan/ && $3+0==$3 {sum+=$3; count++} END {print (count>0 ? sum/count : "NA")}' ${OUT_PATH}/${OUT_PREFIX}_polyA.sites.pi)

# Calculate for polyB
# vcftools --gzvcf ${VCF_onlyLGAUTOSOMES} --keep "$polyB_popfile" --site-pi --out ${OUT_PATH}/${OUT_PREFIX}_polyB
PI_POLYB=$(awk 'NR>1 && $3 !~ /-?nan/ && $3+0==$3 {sum+=$3; count++} END {print (count>0 ? sum/count : "NA")}' ${OUT_PATH}/${OUT_PREFIX}_polyB.sites.pi)

# Calculate for polyC
# vcftools --gzvcf ${VCF_onlyLGAUTOSOMES} --keep "$polyC_popfile" --site-pi --out ${OUT_PATH}/${OUT_PREFIX}_polyC
PI_POLYC=$(awk 'NR>1 && $3 !~ /-?nan/ && $3+0==$3 {sum+=$3; count++} END {print (count>0 ? sum/count : "NA")}' ${OUT_PATH}/${OUT_PREFIX}_polyC.sites.pi)



# Combine results into one file
echo -e "Treatment\tRep\tMean_pi" > ${OUT_PATH}/genomewide_pi_summary.txt
echo -e "GA1\tA\t$PI_GA1" >> ${OUT_PATH}/genomewide_pi_summary.txt
echo -e "mono\tA\t$PI_MONOA" >> ${OUT_PATH}/genomewide_pi_summary.txt
echo -e "mono\tB\t$PI_MONOB" >> ${OUT_PATH}/genomewide_pi_summary.txt
echo -e "mono\tC\t$PI_MONOC" >> ${OUT_PATH}/genomewide_pi_summary.txt
echo -e "poly\tA\t$PI_POLYA" >> ${OUT_PATH}/genomewide_pi_summary.txt
echo -e "poly\tB\t$PI_POLYB" >> ${OUT_PATH}/genomewide_pi_summary.txt
echo -e "poly\tC\t$PI_POLYC" >> ${OUT_PATH}/genomewide_pi_summary.txt

echo "Done! Summary saved to genomewide_pi_summary.txt"