#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 20G
##SBATCH --cpus-per-task 8
#SBATCH -t 150:00:00                               # Set time limit 
#SBATCH --job-name=ss_run.snpEff_on_baypass_outliers       #Job name
#SBATCH -o ss_run.snpEff_on_baypass_outliers.out           #Standard output log
#SBATCH -e ss_run.snpEff_on_baypass_outliers.err           #Standard error log


# # set up the environment
# module add SnpEff/4.3t
module add bcftools
module add htslib/1.15.1
module add SnpEff/4.3t
module add bcftools

MASTER=~/ss
INPUT_VCF=${MASTER}/snpEff/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_GA1_mono_poly__codingRegions.vcf.gz
ANNOTATED_VCF_UNZIPPED=${MASTER}/snpEff/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_GA1_mono_poly__codingRegions_SNPEFF.vcf
ANNOTATED_VCF=${MASTER}/snpEff/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_GA1_mono_poly__codingRegions_SNPEFF.vcf.gz



### locations files
REGIONS_BED="${MASTER}/baypass/R/outlier_output/outlier_BEDs/ss_MonoPoly_autosomes_sig_c2_regions_allSNPS_run.bed"
SNPS_BED="${MASTER}/baypass/R/outlier_output/outlier_BEDs/ss_MonoPoly_autosomes_sig_c2_snps_allSNPS_run.bed"
TOP_SNP_REGION="${MASTER}/baypass/R/outlier_output/outlier_BEDs/ss_MonoPoly_autosomes_sig_c2_top_snp_per_region_allSNPS_run.bed"

### names VCFs subsetted to each method of counting sig SNPs
ANNOTATED_REGIONS_VCF="${MASTER}/snpEff/C2_outliers/SNPs_in_sig_C2_regions_snpeff.vcf.gz"
ANNOTATED_SNPS_VCF="${MASTER}/snpEff/C2_outliers/SNPs_with_sig_C2_snpeff.vcf.gz"
ANNOTATED_TOP_SNPS_PER_REGION_VCF="${MASTER}/snpEff/C2_outliers/SNPs_top_in_each_sig_c2_region_snpeff.vcf.gz"


### names VCFs subsetted to each method of counting sig SNPs, each var cat

ANNOTATED_FULL_VCF_SYN="${MASTER}/snpEff/C2_outliers/SNPs_FULL__snpeffSYNONYMOUS.vcf.gz"
ANNOTATED_FULL_VCF_MIS="${MASTER}/snpEff/C2_outliers/SNPs_FULL__snpeffMISSENSE.vcf.gz"
ANNOTATED_FULL_VCF_NON="${MASTER}/snpEff/C2_outliers/SNPs_FULL__snpeffNONSENSE.vcf.gz"

ANNOTATED_REGIONS_VCF_SYN="${MASTER}/snpEff/C2_outliers/SNPs_in_sig_C2_regions__snpeffSYNONYMOUS.vcf.gz"
ANNOTATED_REGIONS_VCF_MIS="${MASTER}/snpEff/C2_outliers/SNPs_in_sig_C2_regions__snpeffMISSENSE.vcf.gz"
ANNOTATED_REGIONS_VCF_NON="${MASTER}/snpEff/C2_outliers/SNPs_in_sig_C2_regions__snpeffNONSENSE.vcf.gz"

ANNOTATED_SNPS_VCF_SYN="${MASTER}/snpEff/C2_outliers/SNPs_with_sig_C2__snpeffSYNONYMOUS.vcf.gz"
ANNOTATED_SNPS_VCF_MIS="${MASTER}/snpEff/C2_outliers/SNPs_with_sig_C2__snpeffMISSENSE.vcf.gz"
ANNOTATED_SNPS_VCF_NON="${MASTER}/snpEff/C2_outliers/SNPs_with_sig_C2__snpeffNONSENSE.vcf.gz"

ANNOTATED_TOP_SNPS_PER_REGION_VCF_SYN="${MASTER}/snpEff/C2_outliers/SNPs_top_in_each_sig_c2_region__snpeffSYNONYMOUS.vcf.gz"
ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS="${MASTER}/snpEff/C2_outliers/SNPs_top_in_each_sig_c2_region__snpeffMISSENSE.vcf.gz"
ANNOTATED_TOP_SNPS_PER_REGION_VCF_NON="${MASTER}/snpEff/C2_outliers/SNPs_top_in_each_sig_c2_region__snpeffNONSENSE.vcf.gz"



# ### take the annotated, coding regions VCF and pull out the SNPs in locations provided in .bed files
# bcftools view -R ${REGIONS_BED} ${ANNOTATED_VCF}  -Oz -o ${ANNOTATED_REGIONS_VCF}

# bcftools view -R ${SNPS_BED} ${ANNOTATED_VCF}  -Oz -o ${ANNOTATED_SNPS_VCF}

# bcftools view -R ${TOP_SNP_REGION} ${ANNOTATED_VCF}  -Oz -o ${ANNOTATED_TOP_SNPS_PER_REGION_VCF}




# #############################################------------------------------------- ANNOTATE THE VCF WITH snpEff

# # USETHE EXISTING T.CASTANEUM DATABASE WITHIN SNPEFF
# # index the output
# snpEff Tribolium_castaneum ${INPUT_VCF} > ${ANNOTATED_VCF_UNZIPPED}
# bgzip ${ANNOTATED_VCF_UNZIPPED}
# bcftools index ${ANNOTATED_VCF}


# ############################################-------------------------------------  FILTER THE OUTPUT FOR SNP IMPACT


# # ## filter the VCFs for each type of C2 match to sites annotated as each category

### ALL DATA FOR COMPARISON
bcftools view -i 'INFO/ANN[*] ~ "synonymous_variant"' -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${ANNOTATED_VCF} -Oz -o ${ANNOTATED_FULL_VCF_SYN}  ## CHECK THAT I CAN GREP THE ANNOTATED VCF FOR THESE REGEXs


bcftools view -i 'INFO/ANN[*] ~ "missense_variant"' -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${ANNOTATED_VCF} -Oz -o ${ANNOTATED_FULL_VCF_MIS}


bcftools view -i 'INFO/ANN[*] ~ "stop_gained"' -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${ANNOTATED_VCF} -Oz -o ${ANNOTATED_FULL_VCF_NON}



# ### REGIONS
# bcftools view -i 'INFO/ANN[*] ~ "synonymous_variant"' ${ANNOTATED_REGIONS_VCF} -Oz -o ${ANNOTATED_REGIONS_VCF_SYN}  ## CHECK THAT I CAN GREP THE ANNOTATED VCF FOR THESE REGEXs


# bcftools view -i 'INFO/ANN[*] ~ "missense_variant"' ${ANNOTATED_REGIONS_VCF} -Oz -o ${ANNOTATED_REGIONS_VCF_MIS}


# bcftools view -i 'INFO/ANN[*] ~ "stop_gained"' ${ANNOTATED_REGIONS_VCF} -Oz -o ${ANNOTATED_REGIONS_VCF_NON}



# ### SNPS
# bcftools view -i 'INFO/ANN[*] ~ "synonymous_variant"' ${ANNOTATED_SNPS_VCF} -Oz -o ${ANNOTATED_SNPS_VCF_SYN}  ## CHECK THAT I CAN GREP THE ANNOTATED VCF FOR THESE REGEXs


# bcftools view -i 'INFO/ANN[*] ~ "missense_variant"' ${ANNOTATED_SNPS_VCF} -Oz -o ${ANNOTATED_SNPS_VCF_MIS}


# bcftools view -i 'INFO/ANN[*] ~ "stop_gained"' ${ANNOTATED_SNPS_VCF} -Oz -o ${ANNOTATED_SNPS_VCF_NON}



### TOP SNPS per REGION
bcftools view -i 'INFO/ANN[*] ~ "synonymous_variant"' ${ANNOTATED_TOP_SNPS_PER_REGION_VCF} -Oz -o ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_SYN}  ## CHECK THAT I CAN GREP THE ANNOTATED VCF FOR THESE REGEXs


bcftools view -i 'INFO/ANN[*] ~ "missense_variant"' ${ANNOTATED_TOP_SNPS_PER_REGION_VCF} -Oz -o ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS}


bcftools view -i 'INFO/ANN[*] ~ "stop_gained"' ${ANNOTATED_TOP_SNPS_PER_REGION_VCF} -Oz -o ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_NON}





########
########  count the variants in each category to make a contingency table for Exact test

# List VCF files
vcfs="${ANNOTATED_FULL_VCF_SYN} ${ANNOTATED_FULL_VCF_NON} ${ANNOTATED_REGIONS_VCF_SYN} ${ANNOTATED_REGIONS_VCF_NON} ${ANNOTATED_SNPS_VCF_SYN} ${ANNOTATED_SNPS_VCF_NON} ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_SYN} ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_NON}"

# Create a header for table
echo -e "filename\tnum_variants" > vcf_row_counts.tsv

# Loop through files and count variants (excluding header lines)
for f in $vcfs; do
    n=$(bcftools view -H "$f" | wc -l)
    echo -e "$f\t$n" > ${MASTER}/snpEff/C2_outliers/vcf_row_counts.tsv
done



echo "End of script"


