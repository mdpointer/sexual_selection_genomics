library(tidyverse)


args <- commandArgs(trailingOnly = TRUE)

## pass in the arguments from the submission script
popX_name <- args[1]
X_coding_alleleFs <- args[2] # "/Users/mdp/Downloads/polyall_SYNONYMOUS_autosomes_minor_alleleFs.tsv"
X_intergenic_alleleFs <- args[3] #"/Users/mdp/Downloads/polyall_INTERGENIC_autosomes_minor_alleleFs.tsv"
popY_name <- args[4]
Y_coding_alleleFs <-  args[5] #"/Users/mdp/Downloads/monoall_SYNONYMOUS_autosomes_minor_alleleFs.tsv"
Y_intergenic_alleleFs <- args[6] # "/Users/mdp/Downloads/monoall_INTERGENIC_autosomes_minor_alleleFs.tsv"
var_category <- args[7]
chrom <- args[8]



    
##################  ARRANGE THE DATA
##################
### combine single comparisons (popX popY variant_category) into single dataframe - CODING
dd_coding_allSNP <- full_join(
    read_tsv(X_coding_alleleFs, col_names=F) %>% 
      select(X1, X2, X12) %>%
      rename(chrom=X1,
             pos=X2,
             alleleF_x=X12),
    read_tsv(Y_coding_alleleFs, col_names=F) %>% 
      select(X1, X2, X12) %>% 
      rename(chrom=X1,
             pos=X2,
             alleleF_y=X12)
  ) %>%
  mutate( chrom = substr(chrom, 3, nchar(chrom))) %>% 
  mutate( chrom = case_when(chrom != 10 ~ paste(0, chrom, sep=""),
                            TRUE ~ chrom)) %>% 
  mutate(alleleF_x = as.numeric(alleleF_x),
       alleleF_y = as.numeric(alleleF_y)) %>% 
  drop_na(alleleF_x, alleleF_y) %>% 
  filter(!(alleleF_x==0 & alleleF_y==0)) %>% # filter out SNPs not segregating in either pop
  filter(!(alleleF_x==1 & alleleF_y==1)) %>% # filter out SNPs not segregating in either pop
  arrange(chrom, pos)

# new DF including only SNPs that are shared (segragating in BOTH pops)
dd_coding_sharedSNP <- dd_coding_allSNP %>% 
  filter(!(alleleF_x==0 | alleleF_y==0))

### combine single comparisons (popX popY variant_category) into single dataframe - INTERGENIC
dd_intergenic_allSNP <- full_join(
  read_tsv(X_intergenic_alleleFs, col_names=F) %>% 
    select(X1, X2, X12) %>%
    rename(chrom=X1,
           pos=X2,
           alleleF_x=X12),
  read_tsv(Y_intergenic_alleleFs, col_names=F) %>% 
    select(X1, X2, X12) %>% 
    rename(chrom=X1,
           pos=X2,
           alleleF_y=X12)
) %>%
  mutate( chrom = substr(chrom, 3, nchar(chrom))) %>% 
  mutate( chrom = case_when(chrom != 10 ~ paste(0, chrom, sep=""),
                            TRUE ~ chrom)) %>% 
  mutate(alleleF_x = as.numeric(alleleF_x),
         alleleF_y = as.numeric(alleleF_y)) %>% 
  drop_na(alleleF_x, alleleF_y) %>% 
  filter(!(alleleF_x==0 & alleleF_y==0)) %>% # filter out SNPs not segregating in either pop
  filter(!(alleleF_x==0 & alleleF_y==0)) %>% # filter out SNPs not segregating in either pop
  arrange(chrom, pos)

# new DF including only SNPs that are shared (segragating in BOTH pops)
dd_intergenic_sharedSNP <- dd_intergenic_allSNP %>% 
  filter(!(alleleF_x==0 | alleleF_y==0))



# Function to compute Rxy from coding and intergenic datasets
compute_Rxy <- function(coding_df, intergenic_df) {
  x_y_coding <- sum(coding_df$alleleF_x * (1 - coding_df$alleleF_y))
  y_x_coding <- sum(coding_df$alleleF_y * (1 - coding_df$alleleF_x))
  
  x_y_intergenic <- sum(intergenic_df$alleleF_x * (1 - intergenic_df$alleleF_y))
  y_x_intergenic <- sum(intergenic_df$alleleF_y * (1 - intergenic_df$alleleF_x))
  
  (x_y_coding / x_y_intergenic) / (y_x_coding / y_x_intergenic)
}






########  COMPUTE JACKKNIFED RXY FOR ALLSNP DATA - include SNPs segregating in at least 1 of the 2 pops
########

# Sample intergenic SNPs to match coding SNP count once (jackknife excludes rows systematically later)
intergenic_sample_allSNP <- dd_intergenic_allSNP %>% sample_n(nrow(dd_coding_allSNP))
# 
# # JACKKNIFE: systematically exclude each SNP once
# n_jackknife <- nrow(dd_coding_allSNP)
# rxy_jackknife_allSNP <- numeric(n_jackknife)
# 
# for (i in seq_len(n_jackknife)) {
#   # Exclude ith SNP from coding
#   dd_coding_subset <- dd_coding_allSNP[-i, ]
#   
#   # Exclude ith SNP from intergenic sample (same index)
#   dd_intergenic_subset <- intergenic_sample_allSNP[-i, ]
#   
#   # Compute Rxy for this jackknife replicate
#   rxy_jackknife_allSNP[i] <- compute_Rxy(dd_coding_subset, dd_intergenic_subset)
# }
# 
# rxy_jackknife_allSNP_result <- tibble( rxy = rxy_jackknife_allSNP,
#                                 n_jackknife = n_jackknife,
#                                 data_included = "snps segregating in at least 1 pop",
#                                 popX = popX_name,
#                                 popY = popY_name,
#                                 variant_category = var_category,
#                                 chrom = chrom)


# Set the proportion to drop per replicate
jackknife_prop <- 0.10

# Total number of SNPs
n_total <- nrow(dd_coding_allSNP)

# Number of SNPs to drop per replicate
drop_n <- ceiling(n_total * jackknife_prop)

# Number of jackknife replicates
n_jackknife <- floor(n_total / drop_n)

# Preallocate vector to store Rxy values
rxy_jackknife_allSNP <- numeric(n_jackknife)

# Sample indices once, and split into blocks of ~10%
jackknife_indices <- split(sample(1:n_total), ceiling(seq_along(1:n_total) / drop_n))[1:n_jackknife]

# Loop through and calculate Rxy with 10% of SNPs dropped
for (i in seq_len(n_jackknife)) {
  drop_idx <- jackknife_indices[[i]]
  
  dd_coding_subset <- dd_coding_allSNP[-drop_idx, ]
  dd_intergenic_subset <- intergenic_sample_allSNP[-drop_idx, ]
  
  rxy_jackknife_allSNP[i] <- compute_Rxy(dd_coding_subset, dd_intergenic_subset)
}

# Store results
rxy_jackknife_allSNP_result <- tibble(
  rxy = rxy_jackknife_allSNP,
  n_jackknife = n_jackknife,
  prop_dropped = jackknife_prop,
  data_included = "snps segregating in at least 1 pop",
  popX = popX_name,
  popY = popY_name,
  variant_category = var_category,
  chrom = chrom
)




########  COMPUTE JACKKNIFED RXY FOR SHAREDSNP DATA - include only snps segregating in both populations (in coding and intergenic data)
########

# Sample intergenic SNPs to match coding SNP count once (jackknife excludes rows systematically later)
intergenic_sample_sharedSNP <- dd_intergenic_sharedSNP %>% sample_n(nrow(dd_coding_sharedSNP))

# # JACKKNIFE: systematically exclude each SNP once
# n_jackknife <- nrow(dd_coding_sharedSNP)
# rxy_jackknife_sharedSNP <- numeric(n_jackknife)
# 
# for (i in seq_len(n_jackknife)) {
#   # Exclude ith SNP from coding
#   dd_coding_subset <- dd_coding_sharedSNP[-i, ]
#   
#   # Exclude ith SNP from intergenic sample (same index)
#   dd_intergenic_subset <- intergenic_sample_sharedSNP[-i, ]
#   
#   # Compute Rxy for this jackknife replicate
#   rxy_jackknife_sharedSNP[i] <- compute_Rxy(dd_coding_subset, dd_intergenic_subset)
# }
# 
# rxy_jackknife_sharedSNP_result <- tibble( rxy = rxy_jackknife_sharedSNP,
#                                        n_jackknife = n_jackknife,
#                                        data_included = "shared snps segregating in both pops",
#                                        popX = popX_name,
#                                        popY = popY_name,
#                                        variant_category = var_category,
#                                        chrom = chrom)





# Total number of SNPs
n_total <- nrow(dd_coding_sharedSNP)

# Number of SNPs to drop per replicate
drop_n <- ceiling(n_total * jackknife_prop)

# Number of jackknife replicates
n_jackknife <- floor(n_total / drop_n)

# Preallocate vector to store Rxy values
rxy_jackknife_sharedSNP <- numeric(n_jackknife)

# Sample indices once, and split into blocks of ~10%
jackknife_indices <- split(sample(1:n_total), ceiling(seq_along(1:n_total) / drop_n))[1:n_jackknife]

# Loop through and calculate Rxy with 10% of SNPs dropped
for (i in seq_len(n_jackknife)) {
  drop_idx <- jackknife_indices[[i]]
  
  dd_coding_subset <- dd_coding_sharedSNP[-drop_idx, ]
  dd_intergenic_subset <- intergenic_sample_sharedSNP[-drop_idx, ]
  
  rxy_jackknife_sharedSNP[i] <- compute_Rxy(dd_coding_subset, dd_intergenic_subset)
}

# Store results
rxy_jackknife_sharedSNP_result <- tibble(
  rxy = rxy_jackknife_sharedSNP,
  n_jackknife = n_jackknife,
  prop_dropped = jackknife_prop,
  data_included = "snps segregating in at least 1 pop",
  popX = popX_name,
  popY = popY_name,
  variant_category = var_category,
  chrom = chrom
)






out_data <- rbind(
  rxy_jackknife_allSNP_result,
  rxy_jackknife_sharedSNP_result
)


output_filename <- paste("~/ss/Rxy/", popX_name, popY_name, var_category, chrom, ".tsv", sep="_")

out_data %>% 
  write_delim(path.expand(output_filename), delim = "\t", col_names = FALSE)
