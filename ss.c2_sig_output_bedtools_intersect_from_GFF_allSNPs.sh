#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH --mem 20G
##SBATCH --cpus-per-task 8
#SBATCH -t 150:00:00                               # Set time limit 
#SBATCH --job-name=ss.c2_sig_output_bedtools_intersect_from_GFF_allSNPs        #Job name
#SBATCH -o ss.c2_sig_output_bedtools_intersect_from_GFF_allSNPs.out     #Standard output log
#SBATCH -e ss.c2_sig_output_bedtools_intersect_from_GFF_allSNPs.err     #Standard error log


module load bedtools/2.29.2


MASTER=~/ss


# Define the folder containing your SNP files
REGIONS_FOLDER=${MASTER}/baypass/R/outlier_output/outlier_BEDs
GENE_FILE=${MASTER}/data/Tribolium_castaneum.Tcas5.2.59.gff3
OUTPUT_FOLDER=${MASTER}/baypass/R/outlier_output/bedtools_intersect_output

# Sort the GENE_FILE and save it as a temporary sorted file
SORTED_GENE_FILE="${GENE_FILE}.sorted"
bedtools sort -i "$GENE_FILE" > "$SORTED_GENE_FILE"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_FOLDER

# Loop over each SNP file in the folder
for REGIONS_FILE in $REGIONS_FOLDER/*.bed; do
    # Extract the filename without the path or extension
    BASENAME=$(basename $REGIONS_FILE .bed)
    
    # Define the output file path
    OUTPUT_FILE="${OUTPUT_FOLDER}/${BASENAME}_genes.bed"

    # sort the input regions file
    SORTED_REGIONS_FILE="${REGIONS_FILE}.sorted"
    bedtools sort -i "$REGIONS_FILE" > "$SORTED_REGIONS_FILE"
    
    # Run bedtools closest to find the nearest gene for each SNP
    bedtools intersect -a ${SORTED_REGIONS_FILE} -b ${SORTED_GENE_FILE} -wa -wb > ${OUTPUT_FILE}_genesContainingSNPs

    bedtools closest -a ${SORTED_REGIONS_FILE}  -b ${SORTED_GENE_FILE} -d > ${OUTPUT_FILE}_genesClosestToSNPs

    awk '$6 == "gene"' ${OUTPUT_FILE}_genesContainingSNPs > ${OUTPUT_FILE}_genesContainingSNPs_filtered_to_genes
    awk '$6 == "gene"' ${OUTPUT_FILE}_genesClosestToSNPs > ${OUTPUT_FILE}_genesClosestToSNPs_filtered_to_genes

    awk -F'\t' '{split($12, a, ";"); print a[1]}' ${OUTPUT_FILE}_genesContainingSNPs_filtered_to_genes | grep -o 'TC[0-9]*' > ${OUTPUT_FILE}_genesContainingSNPs_geneIDs_only
    awk -F'\t' '{split($12, a, ";"); print a[1]}' ${OUTPUT_FILE}_genesClosestToSNPs_filtered_to_genes | grep -o 'TC[0-9]*' > ${OUTPUT_FILE}_genesClosestToSNPs_geneIDs_only

    rm ${SORTED_REGIONS_FILE}
done

echo "Bedtools -intersect....DONE"



#########
bedtools closest -a ~/ss/baypass/R/outlier_output/outlier_BEDs/ss_MonoPoly_autosomes_sig_c2_snps_allSNPS_run_sorted.bed -b ~/ss/data/Tribolium_castaneum.Tcas5.2.59.gff3.sorted -d > TEST_snps_genesClosestToSNPs

bedtools closest -a ~/ss/baypass/R/outlier_output/outlier_BEDs/ss_MonoPoly_autosomes_sig_c2_top_snp_per_region_allSNPS_run_sorted.bed -b ~/ss/data/Tribolium_castaneum.Tcas5.2.59.gff3.sorted_genesOnly -d > TEST_TOPsnps_genesClosestToSNPs



