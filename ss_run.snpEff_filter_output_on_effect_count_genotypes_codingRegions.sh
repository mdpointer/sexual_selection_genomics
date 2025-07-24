#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 20G
##SBATCH --cpus-per-task 8
#SBATCH -t 150:00:00                               # Set time limit 
#SBATCH --job-name=ss_run.snpEff_codingRegions_filter_output_on_effect_count_genotypes       #Job name
#SBATCH -o ss_run.snpEff_codingRegions_filter_output_on_effect_count_genotypes.out           #Standard output log
#SBATCH -e ss_run.snpEff_codingRegions_filter_output_on_effect_count_genotypes.err           #Standard error log


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

HIGH_IMPACT_SNP_VCF=${MASTER}/snpEff/codingRegions_snpEff_HIGHimpact.vcf.gz
MODERATE_IMPACT_SNP_VCF=${MASTER}/snpEff/codingRegions_snpEff_MODERATEimpact.vcf.gz
LOW_IMPACT_SNP_VCF=${MASTER}/snpEff/codingRegions_snpEff_LOWimpact.vcf.gz
MODIFIER_IMPACT_SNP_VCF=${MASTER}/snpEff/codingRegions_snpEff_MODIFIERimpact.vcf.gz
SYNONYMOUS_SNP_VCF=${MASTER}/snpEff/codingRegions_snpEff_SYNONYMOUS.vcf.gz
MISSENSE_SNP_VCF=${MASTER}/snpEff/codingRegions_snpEff_MISSENSE.vcf.gz
NONSENSE_SNP_VCF=${MASTER}/snpEff/codingRegions_snpEff_NONSENSE.vcf.gz

SAMPLE_NAMES_FILE=${MASTER}/snpEff/codingRegions_sample_names_file.txt

ALL_SNP_STATS_OUTPUT_AUTOSOMES=${MASTER}/snpEff/codingRegions_all_snp_stats_autosomes.txt
ALL_SNP_STATS_OUTPUT_X=${MASTER}/snpEff/codingRegions_all_snp_stats_X.txt

HIGH_IMPACT_SNP_STATS_OUTPUT_AUTOSOMES=${MASTER}/snpEff/codingRegions_high_impact_stats_autosomes.txt
HIGH_IMPACT_SNP_STATS_OUTPUT_X=${MASTER}/snpEff/codingRegions_high_impact_stats_X.txt

MODERATE_IMPACT_SNP_STATS_OUTPUT_AUTOSOMES=${MASTER}/snpEff/codingRegions_moderate_impact_stats_autosomes.txt
MODERATE_IMPACT_SNP_STATS_OUTPUT_X=${MASTER}/snpEff/codingRegions_moderate_impact_stats_X.txt

LOW_IMPACT_SNP_STATS_OUTPUT_AUTOSOMES=${MASTER}/snpEff/codingRegions_low_impact_stats_autosomes.txt
LOW_IMPACT_SNP_STATS_OUTPUT_X=${MASTER}/snpEff/codingRegions_low_impact_stats_X.txt

MODIFIER_IMPACT_SNP_STATS_OUTPUT_AUTOSOMES=${MASTER}/snpEff/codingRegions_modifier_impact_stats_autosomes.txt
MODIFIER_IMPACT_SNP_STATS_OUTPUT_X=${MASTER}/snpEff/codingRegions_modifier_impact_stats_X.txt

SYNONYMOUS_SNP_STATS_OUTPUT_AUTOSOMES=${MASTER}/snpEff/codingRegions_synonymous_snp_stats_autosomes.txt
SYNONYMOUS_SNP_STATS_OUTPUT_X=${MASTER}/snpEff/codingRegions_synonymous_snp_stats_X.txt

MISSENSE_SNP_STATS_OUTPUT_AUTOSOMES=${MASTER}/snpEff/codingRegions_missense_snp_stats_autosomes.txt
MISSENSE_SNP_STATS_OUTPUT_X=${MASTER}/snpEff/codingRegions_missense_snp_stats_X.txt

NONSENSE_SNP_STATS_OUTPUT_AUTOSOMES=${MASTER}/snpEff/codingRegions_nonsense_snp_stats_autosomes.txt
NONSENSE_SNP_STATS_OUTPUT_X=${MASTER}/snpEff/codingRegions_nonsense_snp_stats_X.txt



# #############################################------------------------------------- ANNOTATE THE VCF WITH snpEff

# # USETHE EXISTING T.CASTANEUM DATABASE WITHIN SNPEFF
# # index the output
# snpEff Tribolium_castaneum ${INPUT_VCF} > ${ANNOTATED_VCF_UNZIPPED}
# bgzip ${ANNOTATED_VCF_UNZIPPED}
# bcftools index ${ANNOTATED_VCF}


# ############################################-------------------------------------  FILTER THE OUTPUT FOR SNP IMPACT


# # ## filter the VCF to sites annotated as HIGH or MODERATE|HIGH
# # ## index the output
# bcftools view -i 'INFO/ANN[*] ~ "HIGH"' ${ANNOTATED_VCF} -Oz -o ${HIGH_IMPACT_SNP_VCF}
# bcftools index ${HIGH_IMPACT_SNP_VCF}

# bcftools view -i 'INFO/ANN[*] ~ "MODERATE"' ${ANNOTATED_VCF} -Oz -o ${MODERATE_IMPACT_SNP_VCF}
# bcftools index ${MODERATE_IMPACT_SNP_VCF}

# bcftools view -i 'INFO/ANN[*] ~ "LOW"' ${ANNOTATED_VCF} -Oz -o ${LOW_IMPACT_SNP_VCF}
# bcftools index ${LOW_IMPACT_SNP_VCF}

# bcftools view -i 'INFO/ANN[*] ~ "MODIFIER"' ${ANNOTATED_VCF} -Oz -o ${MODIFIER_IMPACT_SNP_VCF}
# bcftools index ${MODIFIER_IMPACT_SNP_VCF}

# bcftools view -i 'INFO/ANN[*] ~ "synonymous_variant"' ${ANNOTATED_VCF} -Oz -o ${SYNONYMOUS_SNP_VCF}  ## CHECK THAT I CAN GREP THE ANNOTATED VCF FOR THESE REGEXs
# bcftools index ${SYNONYMOUS_SNP_VCF}

# bcftools view -i 'INFO/ANN[*] ~ "missense_variant"' ${ANNOTATED_VCF} -Oz -o ${MISSENSE_SNP_VCF}
# bcftools index ${MISSENSE_SNP_VCF}

# bcftools view -i 'INFO/ANN[*] ~ "stop_gained"' ${ANNOTATED_VCF} -Oz -o ${NONSENSE_SNP_VCF}
# bcftools index ${NONSENSE_SNP_VCF}




# #############################################-------------------------------------  GET SAMPLE-LEVEL GENOTYPE COUNTS FOR RELEVANT VARIANTS

# # get the sample names from the VCF
# bcftools query -l ${ANNOTATED_VCF} > ${SAMPLE_NAMES_FILE}


# ### Get the variant stats for each individual
# # NB: the -A argument in the grep should be the number of samples +1

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${ANNOTATED_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${ALL_SNP_STATS_OUTPUT_AUTOSOMES}
# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${ANNOTATED_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${ALL_SNP_STATS_OUTPUT_X}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${HIGH_IMPACT_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${HIGH_IMPACT_SNP_STATS_OUTPUT_AUTOSOMES}
# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${HIGH_IMPACT_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${HIGH_IMPACT_SNP_STATS_OUTPUT_X}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${MODERATE_IMPACT_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${MODERATE_IMPACT_SNP_STATS_OUTPUT_AUTOSOMES}
# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${MODERATE_IMPACT_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${MODERATE_IMPACT_SNP_STATS_OUTPUT_X}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${LOW_IMPACT_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${LOW_IMPACT_SNP_STATS_OUTPUT_AUTOSOMES}
# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${LOW_IMPACT_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${LOW_IMPACT_SNP_STATS_OUTPUT_X}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${MODIFIER_IMPACT_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${MODIFIER_IMPACT_SNP_STATS_OUTPUT_AUTOSOMES}
# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${MODIFIER_IMPACT_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${MODIFIER_IMPACT_SNP_STATS_OUTPUT_X}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${SYNONYMOUS_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${SYNONYMOUS_SNP_STATS_OUTPUT_AUTOSOMES}
# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${SYNONYMOUS_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${SYNONYMOUS_SNP_STATS_OUTPUT_X}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${MISSENSE_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${MISSENSE_SNP_STATS_OUTPUT_AUTOSOMES}
# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${MISSENSE_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${MISSENSE_SNP_STATS_OUTPUT_X}

# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LG2,LG3,LG4,LG5,LG6,LG7,LG8,LG9,LG10 ${NONSENSE_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${NONSENSE_SNP_STATS_OUTPUT_AUTOSOMES}
# bcftools stats -S ${SAMPLE_NAMES_FILE} -r LGX ${NONSENSE_SNP_VCF} | grep -A 109 "^# PSC, Per-sample counts" > ${NONSENSE_SNP_STATS_OUTPUT_X}







# Extract genotypes from LGX and summarize per sample: handles haploid (0, 1) and diploid (0/0, 0/1, 1/1)
bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${ANNOTATED_VCF} | \
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
}' > ${ALL_SNP_STATS_OUTPUT_X}

bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${HIGH_IMPACT_SNP_VCF} | \
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
}' > ${HIGH_IMPACT_SNP_STATS_OUTPUT_X}


bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${MODERATE_IMPACT_SNP_VCF} | \
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
}' > ${MODERATE_IMPACT_SNP_STATS_OUTPUT_X}

bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${LOW_IMPACT_SNP_VCF} | \
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
}' > ${LOW_IMPACT_SNP_STATS_OUTPUT_X}


bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${MODIFIER_IMPACT_SNP_VCF} | \
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
}' > ${MODIFIER_IMPACT_SNP_STATS_OUTPUT_X}

bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${SYNONYMOUS_SNP_VCF} | \
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
}' > ${SYNONYMOUS_SNP_STATS_OUTPUT_X}

bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${MISSENSE_SNP_VCF} | \
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
}' > ${MISSENSE_SNP_STATS_OUTPUT_X}


bcftools query -r LGX -f '[%SAMPLE=%GT\t]\n' ${NONSENSE_SNP_VCF} | \
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
}' > ${NONSENSE_SNP_STATS_OUTPUT_X}





echo "End of script"