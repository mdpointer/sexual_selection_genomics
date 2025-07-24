#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 10G
##SBATCH --cpus-per-task 4
#SBATCH -t 150:00:00                               # Set time limit 
#SBATCH --job-name=ss_run.SIFT_on_baypass_outliers       #Job name
#SBATCH -o ss_run.SIFT_on_baypass_outliers.out           #Standard output log
#SBATCH -e ss_run.SIFT_on_baypass_outliers.err           #Standard error log


# # set up the environment
module add SIFT4G/2.0.0
module add bcftools
module add htslib/1.15.1

MASTER=~/ss

# begin with VCFs annotated as missense and subsetted to sig baypass outliers in the SnpEff script
ANNOTATED_FULL_VCF_MIS="${MASTER}/snpEff/C2_outliers/SNPs_FULL__snpeffMISSENSE.vcf.gz"
ANNOTATED_FULL_VCF_MIS_UNZIPPED="${MASTER}/snpEff/C2_outliers/SNPs_FULL__snpeffMISSENSE.vcf"
ANNOTATED_REGIONS_VCF_MIS="${MASTER}/snpEff/C2_outliers/SNPs_in_sig_C2_regions__snpeffMISSENSE.vcf.gz"
ANNOTATED_REGIONS_VCF_MIS_UNZIPPED="${MASTER}/snpEff/C2_outliers/SNPs_in_sig_C2_regions__snpeffMISSENSE.vcf"
ANNOTATED_SNPS_VCF_MIS="${MASTER}/snpEff/C2_outliers/SNPs_with_sig_C2__snpeffMISSENSE.vcf.gz"
ANNOTATED_SNPS_VCF_MIS_UNZIPPED="${MASTER}/snpEff/C2_outliers/SNPs_with_sig_C2__snpeffMISSENSE.vcf"
ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS="${MASTER}/snpEff/C2_outliers/SNPs_top_in_each_sig_c2_region__snpeffMISSENSE.vcf.gz"
ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_UNZIPPED="${MASTER}/snpEff/C2_outliers/SNPs_top_in_each_sig_c2_region__snpeffMISSENSE.vcf"

ANNOTATED_FULL_VCF_MIS_SIFT="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_FULL__snpeffMISSENSE_SIFTpredictions.vcf.gz"
ANNOTATED_REGIONS_VCF_MIS_SIFT="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_in_sig_C2_regions__snpeffMISSENSE_SIFTpredictions.vcf.gz"
ANNOTATED_SNPS_VCF_MIS_SIFT="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_with_sig_C2__snpeffMISSENSE_SIFTpredictions.vcf.gz"
ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_top_in_each_sig_c2_region__snpeffMISSENSE_SIFTpredictions.vcf.gz"

ANNOTATED_FULL_VCF_MIS_SIFT_DEL="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_FULL__snpeffMISSENSE_SIFTpredictions_DELETERIOUS.vcf.gz"
ANNOTATED_FULL_VCF_MIS_SIFT_TOL="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_FULL__snpeffMISSENSE_SIFTpredictions_TOLERATED.vcf.gz"
ANNOTATED_FULL_VCF_MIS_SIFT_DEL_HIGHCONF="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_FULL__snpeffMISSENSE_SIFTpredictions_DELETERIOUShighConf.vcf.gz"

ANNOTATED_REGIONS_VCF_MIS_SIFT_DEL="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_in_sig_C2_regions__snpeffMISSENSE_SIFTpredictions_DELETERIOUS.vcf.gz"
ANNOTATED_REGIONS_VCF_MIS_SIFT_TOL="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_in_sig_C2_regions__snpeffMISSENSE_SIFTpredictions_TOLERATED.vcf.gz"
ANNOTATED_REGIONS_VCF_MIS_SIFT_DEL_HIGHCONF="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_in_sig_C2_regions__snpeffMISSENSE_SIFTpredictions_DELETERIOUShighConf.vcf.gz"

ANNOTATED_SNPS_VCF_MIS_SIFT_DEL="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_with_sig_C2__snpeffMISSENSE_SIFTpredictions_DELETERIOUS.vcf.gz"
ANNOTATED_SNPS_VCF_MIS_SIFT_TOL="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_with_sig_C2__snpeffMISSENSE_SIFTpredictions_TOLERATED.vcf.gz"
ANNOTATED_SNPS_VCF_MIS_SIFT_DEL_HIGHCONF="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_with_sig_C2__snpeffMISSENSE_SIFTpredictions_DELETERIOUS_highConf.vcf.gz"

ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT_DEL="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_top_in_each_sig_c2_region__snpeffMISSENSE_SIFTpredictions_DELETERIOUS.vcf.gz"
ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT_TOL="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_top_in_each_sig_c2_region__snpeffMISSENSE_SIFTpredictions_TOLERATED.vcf.gz"
ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT_DEL_HIGHCONF="${MASTER}/SIFT/SIFT_for_baypass_outliers/SNPs_top_in_each_sig_c2_region__snpeffMISSENSE_SIFTpredictions_DELETERIOUShighConf.vcf.gz"





RESULTS_FOLDER=${MASTER}/SIFT/SIFT_for_baypass_outliers


# #############################################------------------------------------- ANNOTATE THE VCFs WITH SIFT - one per set of input sites

# gunzip ${ANNOTATED_FULL_VCF_MIS}

# java -jar ${MASTER}/SIFT/SIFT4G_Annotator.jar -t 4 -c -i ${ANNOTATED_FULL_VCF_MIS_UNZIPPED} -d ${MASTER}/SIFT/SIFT_database_Tcas5.2 -r ${RESULTS_FOLDER} -t

# bgzip ${ANNOTATED_FULL_VCF_MIS_UNZIPPED}
# bcftools index --csi ${ANNOTATED_FULL_VCF_MIS_SIFT}



# gunzip ${ANNOTATED_REGIONS_VCF_MIS}

# java -jar ${MASTER}/SIFT/SIFT4G_Annotator.jar -t 4 -c -i ${ANNOTATED_REGIONS_VCF_MIS_UNZIPPED} -d ${MASTER}/SIFT/SIFT_database_Tcas5.2 -r ${RESULTS_FOLDER} -t

# bgzip ${ANNOTATED_REGIONS_VCF_MIS_UNZIPPED}
# bcftools index --csi ${ANNOTATED_REGIONS_VCF_MIS_SIFT}



# gunzip ${ANNOTATED_SNPS_VCF_MIS}

# java -jar ${MASTER}/SIFT/SIFT4G_Annotator.jar -t 4 -c -i ${ANNOTATED_SNPS_VCF_MIS_UNZIPPED} -d ${MASTER}/SIFT/SIFT_database_Tcas5.2 -r ${RESULTS_FOLDER} -t

# bgzip ${ANNOTATED_SNPS_VCF_MIS_UNZIPPED}
# bcftools index --csi ${ANNOTATED_SNPS_VCF_MIS_SIFT}


###### NB - THIS CATEGORY WAS EMPTY COMING FROM SNPEFF AND MAKES THIS SCRIPT HANG, SO COMMENTED OUT ALONG WITH SUBSEQUENT STEPS FOR THIS SET OF LOCATIONS
## gunzip ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS}
#
## java -jar ${MASTER}/SIFT/SIFT4G_Annotator.jar -t 4 -c -i ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_UNZIPPED} -d ${MASTER}/SIFT/SIFT_database_Tcas5.2 -r ${RESULTS_FOLDER} -t
#
## bgzip ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_UNZIPPED}
# bcftools index --csi ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT}


# # #############################################-------------------------------------  FILTER THE ANNOTATED VCF TO GET THE SUBSET VCFs - for each set of input sites

# # create a new VCFs (from the unfiltered VCF) containing variants annotated as deleterious OR tolerated
# bcftools view -i 'INFO/SIFTINFO[*] ~ "DELETERIOUS"' -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${ANNOTATED_FULL_VCF_MIS_SIFT} -Oz -o ${ANNOTATED_FULL_VCF_MIS_SIFT_DEL}
# bcftools index ${ANNOTATED_FULL_VCF_MIS_SIFT_DEL}

# bcftools view -i 'INFO/SIFTINFO[*] ~ "TOLERATED"' -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${ANNOTATED_FULL_VCF_MIS_SIFT} -Oz -o ${ANNOTATED_FULL_VCF_MIS_SIFT_TOL}
# bcftools index ${ANNOTATED_FULL_VCF_MIS_SIFT_TOL}



# bcftools view -i 'INFO/SIFTINFO[*] ~ "DELETERIOUS"' ${ANNOTATED_REGIONS_VCF_MIS_SIFT} -Oz -o ${ANNOTATED_REGIONS_VCF_MIS_SIFT_DEL}
# bcftools index ${ANNOTATED_REGIONS_VCF_MIS_SIFT_DEL}

# bcftools view -i 'INFO/SIFTINFO[*] ~ "TOLERATED"' ${ANNOTATED_REGIONS_VCF_MIS_SIFT} -Oz -o ${ANNOTATED_REGIONS_VCF_MIS_SIFT_TOL}
# bcftools index ${ANNOTATED_REGIONS_VCF_MIS_SIFT_TOL}



# bcftools view -i 'INFO/SIFTINFO[*] ~ "DELETERIOUS"' ${ANNOTATED_SNPS_VCF_MIS_SIFT} -Oz -o ${ANNOTATED_SNPS_VCF_MIS_SIFT_DEL}
# bcftools index ${ANNOTATED_SNPS_VCF_MIS_SIFT_DEL}

# bcftools view -i 'INFO/SIFTINFO[*] ~ "TOLERATED"' ${ANNOTATED_SNPS_VCF_MIS_SIFT} -Oz -o ${ANNOTATED_SNPS_VCF_MIS_SIFT_TOL}
# bcftools index ${ANNOTATED_SNPS_VCF_MIS_SIFT_TOL}


# ## bcftools view -i 'INFO/SIFTINFO[*] ~ "DELETERIOUS"' ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT} -Oz -o ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT_DEL}
# ## bcftools index ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT_DEL}
# #
# ## bcftools view -i 'INFO/SIFTINFO[*] ~ "TOLERATED"' ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT} -Oz -o ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT_TOL}
# ## bcftools index ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT_TOL}

# # # # #############################################------------------------------------- CREATE VERSIONS OF DELETERIOIUS FILES WITHOUT LOW CONFIDENCE CALLS

# bcftools view -i 'INFO/SIFTINFO[*] !~ "Low confidence"' -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${ANNOTATED_FULL_VCF_MIS_SIFT_DEL} -Oz -o ${ANNOTATED_FULL_VCF_MIS_SIFT_DEL_HIGHCONF}
# bcftools index ${ANNOTATED_FULL_VCF_MIS_SIFT_DEL_HIGHCONF}

# bcftools view -i 'INFO/SIFTINFO[*] !~ "Low confidence"' ${ANNOTATED_REGIONS_VCF_MIS_SIFT_DEL} -Oz -o ${ANNOTATED_REGIONS_VCF_MIS_SIFT_DEL_HIGHCONF}
# bcftools index ${ANNOTATED_REGIONS_VCF_MIS_SIFT_DEL_HIGHCONF}

# bcftools view -i 'INFO/SIFTINFO[*] !~ "Low confidence"' ${ANNOTATED_SNPS_VCF_MIS_SIFT_DEL} -Oz -o ${ANNOTATED_SNPS_VCF_MIS_SIFT_DEL_HIGHCONF}
# bcftools index ${ANNOTATED_SNPS_VCF_MIS_SIFT_DEL_HIGHCONF}

# ##bcftools view -i 'INFO/SIFTINFO[*] !~ "Low confidence"' ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT_DEL} -Oz -o ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT_DEL_HIGHCONF}
# ##bcftools index ${ANNOTATED_TOP_SNPS_PER_REGION_VCF_MIS_SIFT_DEL_HIGHCONF}









# ########
# ########  count the variants in each category to make a contingency table for Exact test

# List VCF files
vcfs="${ANNOTATED_FULL_VCF_MIS_SIFT_DEL} ${ANNOTATED_FULL_VCF_MIS_SIFT_TOL} ${ANNOTATED_FULL_VCF_MIS_SIFT_DEL_HIGHCONF} ${ANNOTATED_REGIONS_VCF_MIS_SIFT_DEL} ${ANNOTATED_REGIONS_VCF_MIS_SIFT_TOL} ${ANNOTATED_REGIONS_VCF_MIS_SIFT_DEL_HIGHCONF} ${ANNOTATED_SNPS_VCF_MIS_SIFT_DEL} ${ANNOTATED_SNPS_VCF_MIS_SIFT_TOL} ${ANNOTATED_SNPS_VCF_MIS_SIFT_DEL_HIGHCONF}"

# Create a header for table
echo -e "filename\tnum_variants" > vcf_row_counts.tsv

# Loop through files and count variants (excluding header lines)
for f in $vcfs; do
    n=$(bcftools view -H "$f" | wc -l)
    echo -e "$f\t$n" >> ${RESULTS_FOLDER}/vcf_row_counts.tsv
done



echo "End of script"
