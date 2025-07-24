#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 10G
##SBATCH --cpus-per-task 4
#SBATCH -t 150:00:00                               # Set time limit 
#SBATCH --job-name=ss_run.SIFT_filter_output_on_effect_count_genotypes_codingRegions       #Job name
#SBATCH -o ss_run.SIFT_filter_output_on_effect_count_genotypes_codingRegions.out           #Standard output log
#SBATCH -e ss_run.SIFT_filter_output_on_effect_count_genotypes_codingRegions.err           #Standard error log


# # set up the environment
module add SIFT4G/2.0.0
module add bcftools
module add htslib/1.15.1

MASTER=~/ss
# INPUT_GZVCF=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_GA1_mono_poly__codingRegions.vcf.gz
INPUT_VCF=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_GA1_mono_poly__codingRegions.vcf
# ANNOTATED_VCF_UNZIPPED=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_GA1_mono_poly__codingRegions_SIFTpredictions.vcf
ANNOTATED_VCF=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_GA1_mono_poly__codingRegions_SIFTpredictions.vcf.gz

MISSENSE_SNP_VCF_FROM_snpEff=${MASTER}/snpEff/codingRegions_snpEff_MISSENSE.vcf.gz
MISSENSE_SNP_VCF=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_snpEff_MISSENSE.vcf.gz
MISSENSE_SNP_POSITIONS_FILE=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_snpEff_MISSENSE_snp_positions.txt

MISSENSE_DELETERIOUS_VCF=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_SIFT_MISSENSE_DELETERIOUS.vcf.gz
MISSENSE_DELETERIOUS_VCF_HIGH_CONFIDENCE=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_SIFT_MISSENSE_DELETERIOUS_high_confidence.vcf.gz
MISSENSE_TOLERATED_VCF=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_SIFT_MISSENSE_TOLERATED.vcf.gz

ALL_SNP_DELETERIOUS_VCF=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_SIFT_ALL_SNP_DELETERIOUS.vcf.gz
ALL_SNP_DELETERIOUS_VCF_HIGH_CONFIDENCE=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_SIFT_ALL_SNP_DELETERIOUS_high_confidence.vcf.gz
ALL_SNP_TOLERATED_VCF=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_SIFT_ALL_SNP_TOLERATED.vcf.gz

SAMPLE_NAMES_FILE=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_VCF_sampleNames.txt


ALL_SNP_DELETERIOUS_STATS_OUTPUT_AUTOSOMES=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_all_snp_deleterious_stats_autosomes.txt
ALL_SNP_DELETERIOUS_STATS_OUTPUT_X=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_all_snp_deleterious_stats_X.txt

ALL_SNP_DELETERIOUS_HIGH_CONF_STATS_OUTPUT_AUTOSOMES=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_all_snp_deleterious_high_conf_stats_autosomes.txt
ALL_SNP_DELETERIOUS_HIGH_CONF_STATS_OUTPUT_X=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_all_snp_deleterious_high_conf_stats_X.txt

ALL_SNP_TOLERATED_STATS_OUTPUT_AUTOSOMES=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_all_snp_tolerated_stats_autosomes.txt
ALL_SNP_TOLERATED_STATS_OUTPUT_X=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_all_snp_tolerated_stats_X.txt

MISSENSE_DELETERIOUS_STATS_OUTPUT_AUTOSOMES=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_missense_deleterious_stats_autosomes.txt
MISSENSE_DELETERIOUS_STATS_OUTPUT_X=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_missense_deleterious_snp_stats_X.txt

MISSENSE_DELETERIOUS_HIGH_CONF_STATS_OUTPUT_AUTOSOMES=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_missense_deleterious_high_conf_stats_autosomes.txt
MISSENSE_DELETERIOUS_HIGH_CONF_STATS_OUTPUT_X=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_missense_deleterious_high_conf_snp_stats_X.txt

MISSENSE_TOLERATED_STATS_OUTPUT_AUTOSOMES=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_missense_tolerated_stats_autosomes.txt
MISSENSE_TOLERATED_STATS_OUTPUT_X=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions/codingRegions_missense_tolerated_snp_stats_X.txt




# #############################################------------------------------------- ANNOTATE THE VCF WITH SIFT

RESULTS_FOLDER=${MASTER}/SIFT/SIFT_monoPolyGa1_codingRegions

# gunzip ${INPUT_GZVCF}

# # SIFT
# # Annotate the VCF with SIFT
# java -jar ${MASTER}/SIFT/SIFT4G_Annotator.jar -t 4 -c -i ${INPUT_VCF} -d ${MASTER}/SIFT/SIFT_database_Tcas5.2 -r ${RESULTS_FOLDER} -t

# bgzip ${INPUT_VCF}
# bgzip ${ANNOTATED_VCF_UNZIPPED}
# bcftools index ${ANNOTATED_VCF}


# # #############################################-------------------------------------  FILTER THE ANNOTATED VCF TO GET THE SUBSET VCFs

# # copy the snpEff missense VCF into SIFT folder
# cp ${MISSENSE_SNP_VCF_FROM_snpEff} ${MISSENSE_SNP_VCF}
# # pull out the positions of missense variants identified by snpEff (the SIFT output doesn't allow separation of missense and nonsense)
# bcftools query -f '%CHROM\t%POS\n' ${MISSENSE_SNP_VCF} > ${MISSENSE_SNP_POSITIONS_FILE}

# # create a new VCFs (from the unfiltered VCF) containing variants annotated as deleterious OR tolerated
# bcftools view -i 'INFO/SIFTINFO[*] ~ "DELETERIOUS"' ${ANNOTATED_VCF} -Oz -o ${ALL_SNP_DELETERIOUS_VCF}
# bcftools index ${ALL_SNP_DELETERIOUS_VCF}

# bcftools view -i 'INFO/SIFTINFO[*] ~ "TOLERATED"' ${ANNOTATED_VCF} -Oz -o ${ALL_SNP_TOLERATED_VCF}
# bcftools index ${ALL_SNP_TOLERATED_VCF}

# # create a new VCF (from the VCF of misense variants) containing missense variants annotated as deleterious OR tolerated
# bcftools view -R ${MISSENSE_SNP_POSITIONS_FILE} -i 'INFO/SIFTINFO[*] ~ "DELETERIOUS"' ${ANNOTATED_VCF} -Oz -o ${MISSENSE_DELETERIOUS_VCF}
# bcftools index ${MISSENSE_DELETERIOUS_VCF}

# bcftools view -R ${MISSENSE_SNP_POSITIONS_FILE} -i 'INFO/SIFTINFO[*] ~ "TOLERATED"' ${ANNOTATED_VCF} -Oz -o ${MISSENSE_TOLERATED_VCF}
# bcftools index ${MISSENSE_TOLERATED_VCF}

# # # #############################################------------------------------------- CREATE VERSIONS OF DELETERIOIUS FILES WITHOUT LOW CONFIDENCE CALLS

# bcftools view -i 'INFO/SIFTINFO[*] !~ "Low confidence"' ${ALL_SNP_DELETERIOUS_VCF} -Oz -o ${ALL_SNP_DELETERIOUS_VCF_HIGH_CONFIDENCE}
# bcftools index ${ALL_SNP_DELETERIOUS_VCF_HIGH_CONFIDENCE}

# bcftools view -i 'INFO/SIFTINFO[*] !~ "Low confidence"' ${MISSENSE_DELETERIOUS_VCF} -Oz -o ${MISSENSE_DELETERIOUS_VCF_HIGH_CONFIDENCE}
# bcftools index ${MISSENSE_DELETERIOUS_VCF_HIGH_CONFIDENCE}


# # #############################################-------------------------------------  GET SAMPLE-LEVEL GENOTYPE COUNTS FOR RELEVANT VARIANTS

# # # # get the sample names from the VCF
# # bcftools query -l ${ANNOTATED_VCF} > ${SAMPLE_NAMES_FILE}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${ALL_SNP_DELETERIOUS_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${ALL_SNP_DELETERIOUS_STATS_OUTPUT_AUTOSOMES}
# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${ALL_SNP_DELETERIOUS_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${ALL_SNP_DELETERIOUS_STATS_OUTPUT_X}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${ALL_SNP_DELETERIOUS_VCF_HIGH_CONFIDENCE} | grep -A 109 "^# PSC, Per-sample counts" > ${ALL_SNP_DELETERIOUS_HIGH_CONF_STATS_OUTPUT_AUTOSOMES}
# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${ALL_SNP_DELETERIOUS_VCF_HIGH_CONFIDENCE} | grep -A 109 "^# PSC, Per-sample counts" > ${ALL_SNP_DELETERIOUS_HIGH_CONF_STATS_OUTPUT_X}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${ALL_SNP_TOLERATED_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${ALL_SNP_TOLERATED_STATS_OUTPUT_AUTOSOMES}
# # bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${ALL_SNP_TOLERATED_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${ALL_SNP_TOLERATED_STATS_OUTPUT_X}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${MISSENSE_DELETERIOUS_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${MISSENSE_DELETERIOUS_STATS_OUTPUT_AUTOSOMES}
# # bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${MISSENSE_DELETERIOUS_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${MISSENSE_DELETERIOUS_STATS_OUTPUT_X}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${MISSENSE_TOLERATED_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${MISSENSE_TOLERATED_STATS_OUTPUT_AUTOSOMES}
# # bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${MISSENSE_TOLERATED_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${MISSENSE_TOLERATED_STATS_OUTPUT_X}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${MISSENSE_DELETERIOUS_VCF_HIGH_CONFIDENCE} | grep -A 109 "^# PSC, Per-sample counts" > ${MISSENSE_DELETERIOUS_HIGH_CONF_STATS_OUTPUT_AUTOSOMES}
# # bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${MISSENSE_DELETERIOUS_VCF_HIGH_CONFIDENCE} | grep -A 109 "^# PSC, Per-sample counts" > ${MISSENSE_DELETERIOUS_HIGH_CONF_STATS_OUTPUT_X}




# Extract genotypes from LGX and summarize per sample: handles haploid (0, 1) and diploid (0/0, 0/1, 1/1)
bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${ALL_SNP_DELETERIOUS_VCF} | \
awk '
{
  for (i = 1; i <= NF; i++) {
    split($i, arr, "=")
    sample = arr[1]
    gt = arr[2]

    if (gt == "0/0" || gt == "0") refhom[sample]++
    else if (gt == "1/1" || gt == "1") nonrefhom[sample]++
    else if (gt == "0/1" || gt == "1/0") het[sample]++
    # Ignore other genotypes: ., ./., 1/2, etc.
  }
}
END {
  print "sample\tnRefHom\tnNonRefHom\tnHet"
  PROCINFO["sorted_in"] = "@ind_str_asc"  # Sort output by sample name
  for (s in refhom) {
    printf "%s\t%d\t%d\t%d\n", s, refhom[s]+0, nonrefhom[s]+0, het[s]+0
  }
}' > ${ALL_SNP_DELETERIOUS_STATS_OUTPUT_X}


bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${ALL_SNP_DELETERIOUS_VCF_HIGH_CONFIDENCE} | \
awk '
{
  for (i = 1; i <= NF; i++) {
    split($i, arr, "=")
    sample = arr[1]
    gt = arr[2]

    if (gt == "0/0" || gt == "0") refhom[sample]++
    else if (gt == "1/1" || gt == "1") nonrefhom[sample]++
    else if (gt == "0/1" || gt == "1/0") het[sample]++
    # Ignore other genotypes: ., ./., 1/2, etc.
  }
}
END {
  print "sample\tnRefHom\tnNonRefHom\tnHet"
  PROCINFO["sorted_in"] = "@ind_str_asc"  # Sort output by sample name
  for (s in refhom) {
    printf "%s\t%d\t%d\t%d\n", s, refhom[s]+0, nonrefhom[s]+0, het[s]+0
  }
}' > ${ALL_SNP_DELETERIOUS_HIGH_CONF_STATS_OUTPUT_X}


bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${ALL_SNP_TOLERATED_VCF} | \
awk '
{
  for (i = 1; i <= NF; i++) {
    split($i, arr, "=")
    sample = arr[1]
    gt = arr[2]

    if (gt == "0/0" || gt == "0") refhom[sample]++
    else if (gt == "1/1" || gt == "1") nonrefhom[sample]++
    else if (gt == "0/1" || gt == "1/0") het[sample]++
    # Ignore other genotypes: ., ./., 1/2, etc.
  }
}
END {
  print "sample\tnRefHom\tnNonRefHom\tnHet"
  PROCINFO["sorted_in"] = "@ind_str_asc"  # Sort output by sample name
  for (s in refhom) {
    printf "%s\t%d\t%d\t%d\n", s, refhom[s]+0, nonrefhom[s]+0, het[s]+0
  }
}' > ${ALL_SNP_TOLERATED_STATS_OUTPUT_X}


bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${MISSENSE_DELETERIOUS_VCF} | \
awk '
{
  for (i = 1; i <= NF; i++) {
    split($i, arr, "=")
    sample = arr[1]
    gt = arr[2]

    if (gt == "0/0" || gt == "0") refhom[sample]++
    else if (gt == "1/1" || gt == "1") nonrefhom[sample]++
    else if (gt == "0/1" || gt == "1/0") het[sample]++
    # Ignore other genotypes: ., ./., 1/2, etc.
  }
}
END {
  print "sample\tnRefHom\tnNonRefHom\tnHet"
  PROCINFO["sorted_in"] = "@ind_str_asc"  # Sort output by sample name
  for (s in refhom) {
    printf "%s\t%d\t%d\t%d\n", s, refhom[s]+0, nonrefhom[s]+0, het[s]+0
  }
}' > ${MISSENSE_DELETERIOUS_STATS_OUTPUT_X}

bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${MISSENSE_TOLERATED_VCF} | \
awk '
{
  for (i = 1; i <= NF; i++) {
    split($i, arr, "=")
    sample = arr[1]
    gt = arr[2]

    if (gt == "0/0" || gt == "0") refhom[sample]++
    else if (gt == "1/1" || gt == "1") nonrefhom[sample]++
    else if (gt == "0/1" || gt == "1/0") het[sample]++
    # Ignore other genotypes: ., ./., 1/2, etc.
  }
}
END {
  print "sample\tnRefHom\tnNonRefHom\tnHet"
  PROCINFO["sorted_in"] = "@ind_str_asc"  # Sort output by sample name
  for (s in refhom) {
    printf "%s\t%d\t%d\t%d\n", s, refhom[s]+0, nonrefhom[s]+0, het[s]+0
  }
}' > ${MISSENSE_TOLERATED_STATS_OUTPUT_X}

bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${MISSENSE_DELETERIOUS_VCF_HIGH_CONFIDENCE} | \
awk '
{
  for (i = 1; i <= NF; i++) {
    split($i, arr, "=")
    sample = arr[1]
    gt = arr[2]

    if (gt == "0/0" || gt == "0") refhom[sample]++
    else if (gt == "1/1" || gt == "1") nonrefhom[sample]++
    else if (gt == "0/1" || gt == "1/0") het[sample]++
    # Ignore other genotypes: ., ./., 1/2, etc.
  }
}
END {
  print "sample\tnRefHom\tnNonRefHom\tnHet"
  PROCINFO["sorted_in"] = "@ind_str_asc"  # Sort output by sample name
  for (s in refhom) {
    printf "%s\t%d\t%d\t%d\n", s, refhom[s]+0, nonrefhom[s]+0, het[s]+0
  }
}' > ${MISSENSE_DELETERIOUS_HIGH_CONF_STATS_OUTPUT_X}

echo "End of script"

