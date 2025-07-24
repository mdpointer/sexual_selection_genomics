library(tidyverse)

C2_snps_file<- "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/baypass/outlier_output/ss_baypass_datarun_6pop_MonoPoly_bothSEX_onlyX_allREG_summary_contrast.outbaypass_c2_outlier_snps_pod0.999.tsv"


data <- read.table(C2_snps_file, header =T, sep = "\t") %>% 
  select(CHROM, POS, C2_std)

data %>% 
  write_delim("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/baypass/outlier_output/ss_MonoPoly_autosomes_sig_c2_snps_allSNPS_run.tsv",
              delim="\t",
              col_names=F,
              quote="none")
  

nrow(data)
data %>% 
  count(CHROM)


### output a .bed file of these outlier C2 positions
data %>% 
  mutate(CHROM=CHROM,
         start=POS-1,
         end=POS) %>% 
  select(CHROM, start, end) %>% 
  write_delim("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/baypass/outlier_output/ss_MonoPoly_X_sig_c2_snps_allSNPS_run.bed",
              delim="\t",
              col_names=F,
              quote="none")




# Sort the data by chromosome and position
data <- data %>%
  arrange(CHROM, POS)

# Parameters
max_distance <- 50000  # Maximum distance between SNPs within a region

# Initialize variables for storing regions
regions <- list()
current_region <- NULL

# Group SNPs into significant regions
for (chrom in unique(data$CHROM)) {
  chrom_data <- data %>% filter(CHROM == chrom)
  positions <- chrom_data$POS
  
  # Initialize a region
  region_start <- positions[1]
  region_snps <- c(region_start)
  
  for (i in 2:length(positions)) {
    if (positions[i] - positions[i - 1] < max_distance) {
      region_snps <- c(region_snps, positions[i])
    } else {
      # Finalize the previous region if it has ≥2 SNPs
      if (length(region_snps) >= 2) {
        regions <- append(regions, list(data.frame(
          CHROM = chrom,
          START = region_start,
          END = max(region_snps),
          NUM_SNPS = length(region_snps)
        )))
      }
      # Start a new region
      region_start <- positions[i]
      region_snps <- c(region_start)
    }
  }
  
  # Handle the last region
  if (length(region_snps) >= 2) {
    regions <- append(regions, list(data.frame(
      CHROM = chrom,
      START = region_start,
      END = max(region_snps),
      NUM_SNPS = length(region_snps)
    )))
  }
}

# Combine all regions into a single data frame
region_df <- do.call(rbind, regions)
region_df %>% 
  count(CHROM)
# Save the result to a file

region_df %>% 
write_delim("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/baypass/outlier_output/ss_MonoPoly_X_baypass_C2_significant_regions_50kb.tsv",
            delim="/t",
            col_names=F,
            quote = "none")


## output the regions as a .bed file
region_df %>% 
  mutate(CHROM=CHROM,
         start=START-1,
         end=END-1) %>% 
  select(CHROM, start, end) %>% 
  write_delim("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/baypass/outlier_output/ss_MonoPoly_X_sig_c2_regions_allSNPS_run.bed",
              delim="\t",
              col_names=F,
              quote="none")





## pull out just the top snp from a region
# For each region, filter the SNPs from `data` that fall in the region and pick the top one
top_snps <- region_df %>%
  rowwise() %>%
  do({
    region <- .
    snps_in_region <- data %>%
      filter(
        CHROM == region$CHROM,
        POS >= region$START,
        POS <= region$END
      )
    # Pick the SNP with the highest C2_std in this region
    snps_in_region %>%
      slice_max(order_by = C2_std, n = 1)
  }) %>%
  ungroup() %>% 
  select(CHROM,POS)

## top snps as .tsv
top_snps %>% 
  select(CHROM, POS) %>% 
write_delim("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/baypass/outlier_output/ss_MonoPoly_X_sig_c2_top_snp_per_region_allSNPS_run.tsv",
            delim="\t",
            col_names=F,
            quote="none")

## top snps as .bed
top_snps %>%
  mutate(
    chrom = CHROM,
    start = POS - 1,   # BED format is 0-based start
    end = POS          # BED format is 0-based, half-open
  ) %>%
  select(chrom, start, end) %>%
  write_delim("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/baypass/outlier_output/ss_MonoPoly_X_sig_c2_top_snp_per_region_allSNPS_run.bed",
              delim="\t",
              col_names=F,
              quote="none")


# View the result
print(filtered_snps)
