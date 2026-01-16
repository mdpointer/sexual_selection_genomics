# sexual_selection_genomics

We use whole genome resequencing from long-term lines of Tribolium castaneum evolved under manipulated levels of sexual selection, either polyandry or enforced monandry, to test hypotheses around the purging of mutation load by sexual selection.

This project uses the Tcas5.2 reference genome and annotation.

These scripts analyse a pre-made SNP VCF.


# Analysis structure

## Genetic structure and diversity

### plink_PCA/ss_snps_plink_ldpruned_50100.1_plink_pca_182samples.sh
Adds MAF tags to VCF, filters on MAF(0.05). prunes SNPs for linkage. Runs PCA as QC of population structure / sequencing batch effects.
### plink_PCA/ss_genomics_plink_pca.R
Loads, wrangles and visualises the output of ss_snps_plink_ldpruned_50100.1_plink_pca_182samples.sh

### vcftools_diversity/vcftools_div_pi_per_pop.sh
Defines populations as groups of samples and computes neucleotide diversity per population. Combines outputs to single file.
### vcftools_diversity/vcftools_diversity.R
Summarises and visualises the output of vcftools_div_pi_per_pop.sh

### plink_het/ss.genomewide_heterozygosity_plink.sh
Filters VCF to autosomes. Converts VCF > . bed format. Computes genomewide heterozygosity with PLINK
### plink_het/ss.plink_het_plot.R
Loads, wrangles, visualises and models the output of ss.genomewide_heterozygosity_plink.sh

### plink_ROH/ss_plink_ROH.sh
Filters VCF to autosomes. Converts VCF > . bed format. Computes Runs of Homozygosity with PLINK
### plink_ROH/ss_plink_ROH.R
Loads, wrangles, visualises and models the output of ss_plink_ROH.sh


## BayPass

### generate_MonoPoly_onlyLGautosomes_nonEXONIC10kb_genofile_OMEGA.sh
Harshly filter autosomal SNPs for representation across samples to obtain a core set to use in the BayPass run to generate the allele frequency matrix (omegafile). Output the BayPass input format genofile.
### ss.baypass.MonoPoly_generate_onlyLGautosomes_allREGION_genofile_data_run_allSNPs.sh
Output the BayPass input format genofile for all SNPs for the main BayPass data run (autosomes)

### ss.omega_runs_visualise_matrix.R
Create a plot to visualise the omega matrix generated in generate_MonoPoly_onlyLGautosomes_nonEXONIC10kb_genofile_OMEGA.sh
### ss.sim_from_omega_runs.R
From omegafile created with generate_MonoPoly_onlyLGautosomes_nonEXONIC10kb_genofile_OMEGA.sh, generate a approximate neutral pseudo-observed dataset

### ss.baypass_ARRAY_run_on_simulated_omega_run_genofiles.sh
Run BayPass on the pseudo-observed dataset generated in ss.sim_from_omega_runs.R
### ss.baypass_ARRAY_data_runs_all_contrast_runs_allSNPs.sh
Run BayPass on empirical data

### ss.ARRAY_plot_from_data_runs_allSNPs.R
Plot from output of ss.baypass_ARRAY_data_runs_all_contrast_runs_allSNPs.sh
### ss.ARRAY_plot_from_data_runs_allSNPs.sh
Slurm submit ss.ARRAY_plot_from_data_runs_allSNPs.R

### ss.define_peak_regions_allSNP.R
Take raw BayPass outout from ... and define peak regions

### ss.c2_sig_output_bedtools_intersect_from_GFF_allSNPs.sh
Take a set of regions of interest (provided in .bed files) and intersect them with the genome annotation with BEDTOOLS

### combined_sex_X/omega_run/generate_MonoPoly_bothSEX_onlyLGX_nonEXONIC10kb_genofile_OMEGA.sh
As for corresponding script above, but for SNPs on the X
### combined_sex_X/omega_run/ss.baypass_onlyLGX_bothSEX_sim_from_omega_runs.R
As for corresponding script above, but for SNPs on the X
### combined_sex_X/POD_run/ss.baypass_ARRAY_bothSEX_onlyLGX_run_on_simulated_omega_run_genofiles.sh
As in ss.baypass_ARRAY_run_on_simulated_omega_run_genofiles.sh, but for SNPs on the X
### combined_sex_X/data_run/all_SNPs//ss.ARRAY_onlyLGX_bothSEX_plot_from_data_runs_allSNPs.R
As for corresponding script above, but for SNPs on the X
### combined_sex_X/data_run/all_SNPs//ss.ARRAY_onlyLGX_plot_bothSEX_from_data_runs.sh
As for corresponding script above, but for SNPs on the X
### combined_sex_X/data_run/all_SNPs//ss.baypass.MonoPoly_generate_onlyLGX_bothSEX_allREGION_genofile_data_run_allSNPs.sh
As for corresponding script above, but for SNPs on the X
### combined_sex_X/data_run/all_SNPs//ss.baypass_ARRAY_onlyLGX_bothSEX_data_runs_all_contrast_runs.sh
As for corresponding script above, but for SNPs on the X
### combined_sex_X/data_run/all_SNPs//ss.onlyX_define_peak_regions_allSNP.R
As for corresponding script above, but for SNPs on the X


## Genetic load

### ss_run.snpEff_filter_output_on_effect_count_genotypes_codingRegions.sh
Annotate VCF with snpEff. Extract subsets of snps of each variant category.
### ss_run.snpEff_filter_output_on_effect_count_genotypes_codingRegions_onBaypassOutliers.sh
As in ss_run.snpEff_filter_output_on_effect_count_genotypes_codingRegions.sh, but for only a VCF of significant loci from BayPass.

### SIFT_annotate.sh
Annotates SnpEff annotated VCF with SIFT (using a premade database available from their repo).
Subsets annotated VCF to functional variant categories.
Extracts sample-level genotype counts for each category (separately for autosomal and X loci)

### generate_ratio_syn_plots_and_model.R
Load, wrangle, viaualise and model the output of SIFT_annotate.sh

### SIFT_annotate_for_baypass.sh
As in SIFT_annotate.sh, but for only a VCF of significant loci from BayPass.


### ss.output_the_minor_allele_in_GA1_per_var_cat_codingRegions.sh
Extract the minor allele at each SNP in the ancestral population (coding regions).
### ss.output_the_minor_allele_in_GA1_per_var_cat_intergenicRegions.sh
Extract the minor allele at each SNP in the ancestral population (intergenic regions).
### ss.ALLELE_Fs_autosomes_import_and_summarise.R
Load output of above. Summarise, and output a .RDS file for local plotting.
### ss.ALLELE_Fs_autosomes_import_and_summarise.R
Load output of ss.ALLELE_Fs_autosomes_import_and_summarise.R. Generate allele freq spectrum plots


### ss_compute_Rxy_standardised_jackknifed.R
Take allele frequencies per treatment and compute Rxy + block jacknife to generate variance estimate.
### ss_Rxy_summarise_across_jackknifes.R
Load, wrangle and summarise the output of ss_compute_Rxy_standardised_jackknifed.R


## Survival

### Lumley_selB_survival_data.csv
Family level survival data from Lumley et al. 2015 reanalysed in this study in combination with mutation load metrics

### Survival_analysis_using_lumley_2015_data_plots.R
Survival analysis script performing the above



