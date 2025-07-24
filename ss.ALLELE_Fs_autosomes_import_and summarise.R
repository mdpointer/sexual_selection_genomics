library(tidyverse)


# Specify the folder path
folder_path <- "~/ss/allele_Fs"


# List all .tsv files in the folder
file_list <- list.files(path = folder_path, pattern = "minor_alleleFs\\.tsv$", full.names = TRUE)

# Read and bind all .tsv files into a single dataframe
dd <- file_list %>%
  lapply(function(file) {
    df <- read.table(file, header = FALSE, sep = "\t")
    df <- df %>%
      filter(V5 != ".") %>%             # Exclude rows where V5 is "."
      mutate(V5 = as.numeric(V12)) %>%         # Convert V5 to numeric
      filter(V12 != ".") %>%             # Exclude rows where V12 is "."
      mutate(V12 = as.numeric(V12))        # Convert V5 to numeric
    return(df)
  }) %>%
  bind_rows() %>%
  rename(CHROM = V1,
         POS = V2,
         ref_allele = V3,
         alt_allele = V4,
         alt_allele_f = V5,
         pop = V6,
         var_cat = V7,
         chrom = V8,
         chrom_pos1 = V9,
         chrom_pos2 = V10,
         minor_allele = V11,
         minor_allele_f = V12) %>% 
  select(-ref_allele, -alt_allele, -alt_allele_f, -chrom_pos1, -chrom_pos2)

# add the bins into which each AF falls
dd_prop <- dd %>%
  mutate(bin = cut(minor_allele_f, breaks = seq(0, 1, by = 0.1),
                                include.lowest = TRUE, labels = c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1")))


# compute the proportions
dd_prop <- dd_prop %>%
  group_by(pop, var_cat, chrom, bin) %>%
  summarise(n_freqs = n()) %>%
  ungroup() %>%
  group_by(pop, var_cat) %>%
  mutate(total= sum(n_freqs)) %>%
  mutate(prop= n_freqs/total)

# save the data for local plotting (no X11 on HPC)
dd_prop %>%
  write_rds(file="~/ss/allele_Fs/processed_data/processed_minor_allele_F_proportions_for_AFS_plot.rds")
