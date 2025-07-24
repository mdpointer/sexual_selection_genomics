#!/bin/bash
##SBATCH --mail-type=ALL                    #Mail events (NONE, BEGIN, END, FAIL, ALL)
##SBATCH --mail-user=kfe08xgu@uea.ac.uk           # Where to send mail - change <username> to your  userid
#SBATCH -p compute-64-512                          #Which queue to use
#SBATCH -t 36:00:00                               # Set time limit to 36 hours
#SBATCH --job-name=ss_snps_plink_ldpruned_50100.1_plink_pca_182samples        #Job name
#SBATCH --mem 16G
#SBATCH -o ss_snps_plink_MAF_ldpruned_50_10_0.1_plink_pca_108samples.out                    #Standard output log
#SBATCH -e ss_snps_plink_MAF_ldpruned_50_10_0.1_plink_pca_108samples.err                     #Standard error log



#set up environment
module add vcftools/0.1.16
module add bcftools
module add plink/1.90

# define variables
MASTER=~/ss
VCF=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_GA1_mono_poly_biallelic.vcf.gz
VCF_MAF=${MASTER}/data/bcftools__182inds_sexual_selection__Tribolium_castaneum.Tcas5.2.dna_sm.toplevel__3OutliersRemoved_3LowCovRemoved__ploidyAware__20231028.SNPs_GA1_mono_poly_biallelic_MAF.vcf.gz




# take VCF, fill tags and & MAF filter
bcftools +fill-tags ${VCF} -Ou -- -t MAF | \
bcftools view -e 'MAF<0.05' -Oz -o ${VCF_MAF}



# linkage prune - create files:
# snps to retain    [].prune.in
# and snps to drop  [].prune.out
plink \
--vcf ${VCF_MAF} \
--allow-extra-chr \
--set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 \
--out ${MASTER}/plink/ss_snps_108samples_MAF_plink_ldpruned_50_10_0.1


# run a pca using only snps that passed the LD pruning above
plink \
--vcf ${VCF_MAF} \
--double-id \
--allow-extra-chr \
--set-missing-var-ids @:# \
--extract ${MASTER}/plink/ss_snps_108samples_MAF_plink_ldpruned_50_10_0.1.prune.in \
--make-bed \
--pca \
--out ${MASTER}/plink/ss_snps_182samples_plink_MAF_ldpruned_50_10_0.1_PCA




echo "End of script"

