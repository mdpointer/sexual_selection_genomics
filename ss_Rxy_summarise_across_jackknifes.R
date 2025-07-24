library(tidyverse)
# library(cowplot)

# 
# args <- commandArgs(trailingOnly = TRUE)
# 
# filename <- args[1]
# 
# ## read in concatenated Rxy output files
# ## summarise across bootstraps, within single comparisons,within data types (allSNP vs sharedSNP)
# dd <- read_tsv(filename, col_names=F) %>%
#   rename(Rxy= X1,
#          total_jacknifes = X2,
#          prop_dropped = X3,
#          SNPs_included = X4,
#          popX = X5,
#          popY = X6,
#          variant_category = X7,
#          chrom = X8) %>%
#   group_by(SNPs_included, popX, popY, variant_category, chrom) %>%
#   summarise(mean_Rxy = mean(Rxy, na.rm=T),
#             sd = sd(Rxy, na.rm=T),
#             n = n()) %>%
#   mutate( se = sd / sqrt(n))



## output a combined summarised Rxy file
output_filename <- "~/ss/Rxy/processed_output/Rxy_combined_summarised_output.tsv"

dd %>%
  write_delim(path.expand(output_filename), delim = "\t")





###### LOCAL PLOTTING
file <- "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/Rxy/Rxy_combined_summarised_output.tsv"
dd <- read_tsv(file, col_names=T) %>%
  mutate(pops = paste(popX, popY, sep="_")) %>%
  filter(variant_category != "synonymous") %>%
  mutate( variant_category=case_when(variant_category=="nonsense" ~ "Nonsense",
                                     variant_category=="missense_deleterious" ~ "Missense deleterious",
                                     variant_category=="missense_tolerated" ~ "Missense tolerated")) %>%
  mutate(variant_category = str_replace_all(variant_category, " ", "\n")) %>% 
  mutate(variant_category = fct_relevel(variant_category, c("Nonsense", "Missense\ndeleterious", "Missense\ntolerated")))





#OUTPLOT_poly_mono_privateSNP_autosomes <- 
  dd %>%
  filter(chrom == "autosomes") %>%
  filter(pops == "polyall_monoall") %>%
  filter(SNPs_included == "snps segregating in at least 1 pop") %>%
  ggplot(aes(x = mean_Rxy, y = variant_category)) +
  geom_boxplot(width = 0.7) +
  geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  theme_bw() +
  labs(
    x = expression(R[xy] ~ "(Polyandry / Monandry)"),
    y = ""
  ) +
  theme(
    text = element_text(size = 13),
    # plot.margin = margin(5.5, 130, 5.5, 5.5)
  )+
  xlim(0.73, 1.05)
ggsave(
  filename = "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/Rxy/Rxy_poly_mono_privateSNP_autosomes.png",
  plot = OUTPLOT_poly_mono_privateSNP_autosomes,
  device = "png",
  units="mm",
  height=80,
  width=170,
  dpi=300
)


OUTPLOT_poly_mono_sharedSNP_autosomes <- dd %>%
  filter(chrom == "autosomes") %>%
  filter(pops == "polyall_monoall") %>%
  filter(SNPs_included == "shared snps segregating in both pops") %>%
  ggplot(aes(x = mean_Rxy, y = variant_category)) +
  geom_boxplot(width = 0.7) +
  geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  theme_bw() +
  labs(
    x = expression(R[xy] ~ "(Polyandry / Monandry)"),
    y = ""
  ) +
  theme(
    text = element_text(size = 13),
    # plot.margin = margin(5.5, 130, 5.5, 5.5)
  )+
  xlim(0.73, 1.05)
ggsave(
  filename = "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/Rxy/Rxy_poly_mono_sharedSNP_autosomes.png",
  plot = OUTPLOT_poly_mono_sharedSNP_autosomes,
  device = "png",
  units="mm",
  height=80,
  width=170,
  dpi=300
)


OUTPLOT_poly_mono_privateSNP_X <- dd %>%
  filter(chrom == "X") %>%
  filter(pops == "polyall_monoall") %>%
  filter(SNPs_included == "snps segregating in at least 1 pop") %>%
  ggplot(aes(x = mean_Rxy, y = variant_category)) +
  geom_boxplot(width = 0.7) +
  geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  theme_bw() +
  labs(
    x = expression(R[xy] ~ "(Polyandry / Monandry)"),
    y = ""
  ) +
  theme(
    text = element_text(size = 13),
    # plot.margin = margin(5.5, 130, 5.5, 5.5)
  )+
  xlim(0.2, 1.2)
ggsave(
  filename = "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/Rxy/Rxy_poly_mono_privateSNP_X.png",
  plot = OUTPLOT_poly_mono_privateSNP_X,
  device = "png",
  units="mm",
  height=80,
  width=170,
  dpi=300
)


OUTPLOT_poly_mono_sharedSNP_X <- dd %>%
  filter(chrom == "X") %>%
  filter(pops == "polyall_monoall") %>%
  filter(SNPs_included == "shared snps segregating in both pops") %>%
  ggplot(aes(x = mean_Rxy, y = variant_category)) +
  geom_boxplot(width = 0.7) +
  geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
  geom_vline(xintercept = 1, linetype = "dotted") +
  theme_bw() +
  labs(
    x = expression(R[xy] ~ "(Polyandry / Monandry)"),
    y = ""
  ) +
  theme(
    text = element_text(size = 13),
    # plot.margin = margin(5.5, 130, 5.5, 5.5)
  )+
  xlim(0.2, 1.2)
ggsave(
  filename = "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/Rxy/Rxy_poly_mono_sharedSNP_X.png",
  plot = OUTPLOT_poly_mono_sharedSNP_X,
  device = "png",
  units="mm",
  height=80,
  width=170,
  dpi=300
)











# 
# 
# ##### AUTOSOMES
# ##### PRIVATE SNPS
# 
# mono_ga1_privateSNP_autosomes <- dd %>%
#   filter(chrom == "autosomes") %>%
#   filter(pops == "monoall_GA1") %>%
#   filter(SNPs_included == "snps segregating in at least 1 pop") %>%
#   ggplot(aes(x = mean_Rxy, y = variant_category)) +
#   geom_boxplot(width = 0.7) +
#   geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
#   geom_vline(xintercept = 1, linetype = "dotted") +
#   theme_bw() +
#   labs(
#     # x = expression(R[xy] ~ "(population 1 / population 2)"),
#     y = ""
#   ) +
#   theme(
#     text = element_text(size = 13),
#     plot.margin = margin(5.5, 130, 5.5, 5.5),
#     axis.title.x = element_blank()  
#   )+ 
#   xlim(0.73, 1.05)
# 
# # Combine plot and horizontal title on the right
# mono_ga1_privateSNP_autosomes_title <- ggdraw() +
#   draw_plot(mono_ga1_privateSNP_autosomes, x = 0, width = 0.9) +  # leave space for title on right
#   draw_label("Monandry\n
#              vs\n
#              Ancestral", x = 0.65, y = 0.5,
#              hjust = 0, vjust = 0.5, size = 12, angle = 0, lineheight = 0.55)
# 
# 
# ################################################################################################
# 
# poly_ga1_privateSNP_autosomes <- dd %>%
#   filter(chrom == "autosomes") %>%
#   filter(pops == "polyall_GA1") %>%
#   filter(SNPs_included == "snps segregating in at least 1 pop") %>%
#   ggplot(aes(x = mean_Rxy, y = variant_category)) +
#   geom_boxplot(width = 0.7) +
#   geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
#   geom_vline(xintercept = 1, linetype = "dotted") +
#   theme_bw() +
#   labs(
#     # x = expression(R[xy] ~ "(population 1 / population 2)"),
#     y = ""
#   ) +
#   theme(
#     text = element_text(size = 13),
#     plot.margin = margin(5.5, 130, 5.5, 5.5),
#     axis.title.x = element_blank()  
#   )+
#   xlim(0.73, 1.05)
# 
# # Combine plot and horizontal title on the right
# poly_ga1_privateSNP_autosomes_title <- ggdraw() +
#   draw_plot(poly_ga1_privateSNP_autosomes, x = 0, width = 0.9) +  # leave space for title on right
#   draw_label("Polyandry\n
#              vs\n
#              Ancestral", x = 0.65, y = 0.5,
#              hjust = 0, vjust = 0.5, size = 12, angle = 0, lineheight = 0.55)
# 
# 
# ################################################################################################
# 
# poly_mono_privateSNP_autosomes <- dd %>%
#   filter(chrom == "autosomes") %>%
#   filter(pops == "polyall_monoall") %>%
#   filter(SNPs_included == "snps segregating in at least 1 pop") %>%
#   ggplot(aes(x = mean_Rxy, y = variant_category)) +
#   geom_boxplot(width = 0.7) +
#   geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
#   geom_vline(xintercept = 1, linetype = "dotted") +
#   theme_bw() +
#   labs(
#     x = expression(R[xy] ~ "(population 1 / population 2)"),
#     y = ""
#   ) +
#   theme(
#     text = element_text(size = 13),
#     plot.margin = margin(5.5, 130, 5.5, 5.5)
#   )+
#   xlim(0.73, 1.05)
# 
# # Combine plot and horizontal title on the right
# poly_mono_privateSNP_autosomes_title <- ggdraw() +
#   draw_plot(poly_mono_privateSNP_autosomes, x = 0, width = 0.9) +  # leave space for title on right
#   draw_label("Polyandry\n
#              vs\n
#              Monandry", x = 0.65, y = 0.5,
#              hjust = 0, vjust = 0.5, size = 12, angle = 0, lineheight = 0.55)
# 
# 
# Rxy_privateSNP_autosomes <- ggarrange(mono_ga1_privateSNP_autosomes_title,
#           poly_ga1_privateSNP_autosomes_title,
#           poly_mono_privateSNP_autosomes_title,
#           ncol=1)
# 
# ggsave(
#   filename = "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/Rxy/Rxy_privateSNP_autosomes.png",
#   plot = Rxy_privateSNP_autosomes,
#   device = "png",
#   units="mm",
#   height=120,
#   width=170,
#   dpi=300
# )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##### AUTOSOMES
# ##### SHARED SNPS
# 
# 
# mono_ga1_sharedSNP_autosomes <- dd %>%
#   filter(chrom == "autosomes") %>%
#   filter(pops == "monoall_GA1") %>%
#   filter(SNPs_included == "shared snps segregating in both pops") %>%
#   ggplot(aes(x = mean_Rxy, y = variant_category)) +
#   geom_boxplot(width = 0.7) +
#   geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
#   geom_vline(xintercept = 1, linetype = "dotted") +
#   theme_bw() +
#   labs(
#     # x = expression(R[xy] ~ "(population 1 / population 2)"),
#     y = ""
#   ) +
#   theme(
#     text = element_text(size = 13),
#     plot.margin = margin(5.5, 130, 5.5, 5.5),
#     axis.title.x = element_blank()  
#   )+ 
#   xlim(0.73, 1.05)
# 
# # Combine plot and horizontal title on the right
# mono_ga1_sharedSNP_autosomes_title <- ggdraw() +
#   draw_plot(mono_ga1_sharedSNP_autosomes, x = 0, width = 0.9) +  # leave space for title on right
#   draw_label("Monandry\n
#              vs\n
#              Ancestral", x = 0.65, y = 0.5,
#              hjust = 0, vjust = 0.5, size = 12, angle = 0, lineheight = 0.55)
# 
# ################################################################################################
# 
# poly_ga1_sharedSNP_autosomes <- dd %>%
#   filter(chrom == "autosomes") %>%
#   filter(pops == "polyall_GA1") %>%
#   filter(SNPs_included == "shared snps segregating in both pops") %>%
#   ggplot(aes(x = mean_Rxy, y = variant_category)) +
#   geom_boxplot(width = 0.7) +
#   geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
#   geom_vline(xintercept = 1, linetype = "dotted") +
#   theme_bw() +
#   labs(
#     # x = expression(R[xy] ~ "(population 1 / population 2)"),
#     y = ""
#   ) +
#   theme(
#     text = element_text(size = 13),
#     plot.margin = margin(5.5, 130, 5.5, 5.5),
#     axis.title.x = element_blank()  
#   )+
#   xlim(0.73, 1.05)
# 
# # Combine plot and horizontal title on the right
# poly_ga1_sharedSNP_autosomes_title <- ggdraw() +
#   draw_plot(poly_ga1_sharedSNP_autosomes, x = 0, width = 0.9) +  # leave space for title on right
#   draw_label("Polyandry\n
#              vs\n
#              Ancestral", x = 0.65, y = 0.5,
#              hjust = 0, vjust = 0.5, size = 12, angle = 0, lineheight = 0.55)
# 
# ################################################################################################
# 
# poly_mono_sharedSNP_autosomes <- dd %>%
#   filter(chrom == "autosomes") %>%
#   filter(pops == "polyall_monoall") %>%
#   filter(SNPs_included == "shared snps segregating in both pops") %>%
#   ggplot(aes(x = mean_Rxy, y = variant_category)) +
#   geom_boxplot(width = 0.7) +
#   geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
#   geom_vline(xintercept = 1, linetype = "dotted") +
#   theme_bw() +
#   labs(
#     x = expression(R[xy] ~ "(population 1 / population 2)"),
#     y = ""
#   ) +
#   theme(
#     text = element_text(size = 13),
#     plot.margin = margin(5.5, 130, 5.5, 5.5)
#   )+
#   xlim(0.73, 1.05)
# 
# # Combine plot and horizontal title on the right
# poly_mono_sharedSNP_autosomes_title <- ggdraw() +
#   draw_plot(poly_mono_sharedSNP_autosomes, x = 0, width = 0.9) +  # leave space for title on right
#   draw_label("Polyandry\n
#              vs\n
#              Monandry", x = 0.65, y = 0.5,
#              hjust = 0, vjust = 0.5, size = 12, angle = 0, lineheight = 0.55)
# 
# Rxy_sharedSNP_autosomes <- ggarrange(mono_ga1_sharedSNP_autosomes_title,
#                                       poly_ga1_sharedSNP_autosomes_title,
#                                       poly_mono_sharedSNP_autosomes_title,
#                                       ncol=1)
# 
# ggsave(
#   filename = "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/Rxy/Rxy_sharedSNP_autosomes.png",
#   plot = Rxy_sharedSNP_autosomes,
#   device = "png",
#   units="mm",
#   height=120,
#   width=170,
#   dpi=300
# )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##### X
# ##### PRIVATE SNPS
# 
# mono_ga1_privateSNP_X <- dd %>%
#   filter(chrom == "X") %>%
#   filter(pops == "monoall_GA1") %>%
#   filter(SNPs_included == "snps segregating in at least 1 pop") %>%
#   ggplot(aes(x = mean_Rxy, y = variant_category)) +
#   geom_boxplot(width = 0.7) +
#   geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
#   geom_vline(xintercept = 1, linetype = "dotted") +
#   theme_bw() +
#   labs(
#     # x = expression(R[xy] ~ "(population 1 / population 2)"),
#     y = ""
#   ) +
#   theme(
#     text = element_text(size = 13),
#     plot.margin = margin(5.5, 130, 5.5, 5.5),
#     axis.title.x = element_blank()  
#   )+
#   xlim(0.2, 1.6)
# 
# # Combine plot and horizontal title on the right
# mono_ga1_privateSNP_X_title <- ggdraw() +
#   draw_plot(mono_ga1_privateSNP_X, x = 0, width = 0.9) +  # leave space for title on right
#   draw_label("Monandry\n
#              vs\n
#              Ancestral", x = 0.65, y = 0.5,
#              hjust = 0, vjust = 0.5, size = 12, angle = 0, lineheight = 0.55)
# 
# ################################################################################################
# 
# poly_ga1_privateSNP_X <- dd %>%
#   filter(chrom == "X") %>%
#   filter(pops == "polyall_GA1") %>%
#   filter(SNPs_included == "snps segregating in at least 1 pop") %>%
#   ggplot(aes(x = mean_Rxy, y = variant_category)) +
#   geom_boxplot(width = 0.7) +
#   geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
#   geom_vline(xintercept = 1, linetype = "dotted") +
#   theme_bw() +
#   labs(
#     # x = expression(R[xy] ~ "(population 1 / population 2)"),
#     y = ""
#   ) +
#   theme(
#     text = element_text(size = 13),
#     plot.margin = margin(5.5, 130, 5.5, 5.5),
#     axis.title.x = element_blank()  
#   )+
#   xlim(0.2, 1.6)
# 
# # Combine plot and horizontal title on the right
# poly_ga1_privateSNP_X_title <- ggdraw() +
#   draw_plot(poly_ga1_privateSNP_X, x = 0, width = 0.9) +  # leave space for title on right
#   draw_label("Polyandry\n
#              vs\n
#              Ancestral", x = 0.65, y = 0.5,
#              hjust = 0, vjust = 0.5, size = 12, angle = 0, lineheight = 0.55)
# 
# ################################################################################################
# 
# poly_mono_privateSNP_X <- dd %>%
#   filter(chrom == "X") %>%
#   filter(pops == "polyall_monoall") %>%
#   filter(SNPs_included == "snps segregating in at least 1 pop") %>%
#   ggplot(aes(x = mean_Rxy, y = variant_category)) +
#   geom_boxplot(width = 0.7) +
#   geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
#   geom_vline(xintercept = 1, linetype = "dotted") +
#   theme_bw() +
#   labs(
#     x = expression(R[xy] ~ "(population 1 / population 2)"),
#     y = ""
#   ) +
#   theme(
#     text = element_text(size = 13),
#     plot.margin = margin(5.5, 130, 5.5, 5.5)
#   )+
#   xlim(0.2, 1.6)
# 
# # Combine plot and horizontal title on the right
# poly_mono_privateSNP_X_title <- ggdraw() +
#   draw_plot(poly_mono_privateSNP_X, x = 0, width = 0.9) +  # leave space for title on right
#   draw_label("Polyandry\n
#              vs\n
#              Monandry", x = 0.65, y = 0.5,
#              hjust = 0, vjust = 0.5, size = 12, angle = 0, lineheight = 0.55)
# 
# Rxy_privateSNP_X <- ggarrange(mono_ga1_privateSNP_X_title,
#                                       poly_ga1_privateSNP_X_title,
#                                       poly_mono_privateSNP_X_title,
#                                       ncol=1)
# 
# ggsave(
#   filename = "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/Rxy/Rxy_privateSNP_X.png",
#   plot = Rxy_privateSNP_X,
#   device = "png",
#   units="mm",
#   height=120,
#   width=170,
#   dpi=300
# )
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ##### X
# ##### SHARED SNPS
# 
# 
# mono_ga1_sharedSNP_X <- dd %>%
#   filter(chrom == "X") %>%
#   filter(pops == "monoall_GA1") %>%
#   filter(SNPs_included == "shared snps segregating in both pops") %>%
#   ggplot(aes(x = mean_Rxy, y = variant_category)) +
#   geom_boxplot(width = 0.7) +
#   geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
#   geom_vline(xintercept = 1, linetype = "dotted") +
#   theme_bw() +
#   labs(
#     # x = expression(R[xy] ~ "(population 1 / population 2)"),
#     y = ""
#   ) +
#   theme(
#     text = element_text(size = 13),
#     plot.margin = margin(5.5, 130, 5.5, 5.5),
#     axis.title.x = element_blank()  
#   )+ 
#   xlim(0.2, 2.5)
# 
# # Combine plot and horizontal title on the right
# mono_ga1_sharedSNP_X_title <- ggdraw() +
#   draw_plot(mono_ga1_sharedSNP_X, x = 0, width = 0.9) +  # leave space for title on right
#   draw_label("Monandry\n
#              vs\n
#              Ancestral", x = 0.65, y = 0.5,
#              hjust = 0, vjust = 0.5, size = 12, angle = 0, lineheight = 0.55)
# 
# ################################################################################################
# 
# poly_ga1_sharedSNP_X <- dd %>%
#   filter(chrom == "X") %>%
#   filter(pops == "polyall_GA1") %>%
#   filter(SNPs_included == "shared snps segregating in both pops") %>%
#   ggplot(aes(x = mean_Rxy, y = variant_category)) +
#   geom_boxplot(width = 0.7) +
#   geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
#   geom_vline(xintercept = 1, linetype = "dotted") +
#   theme_bw() +
#   labs(
#     # x = expression(R[xy] ~ "(population 1 / population 2)"),
#     y = ""
#   ) +
#   theme(
#     text = element_text(size = 13),
#     plot.margin = margin(5.5, 130, 5.5, 5.5),
#     axis.title.x = element_blank()  
#   )+
#   xlim(0.2, 2.5)
# 
# # Combine plot and horizontal title on the right
# poly_ga1_sharedSNP_X_title <- ggdraw() +
#   draw_plot(poly_ga1_sharedSNP_X, x = 0, width = 0.9) +  # leave space for title on right
#   draw_label("Polyandry\n
#              vs\n
#              Ancestral", x = 0.65, y = 0.5,
#              hjust = 0, vjust = 0.5, size = 12, angle = 0, lineheight = 0.55)
# 
# ################################################################################################
# 
# poly_mono_sharedSNP_X <- dd %>%
#   filter(chrom == "X") %>%
#   filter(pops == "polyall_monoall") %>%
#   filter(SNPs_included == "shared snps segregating in both pops") %>%
#   ggplot(aes(x = mean_Rxy, y = variant_category)) +
#   geom_boxplot(width = 0.7) +
#   geom_errorbarh(aes(xmin = mean_Rxy - se, xmax = mean_Rxy + se), height = 0.2, color = "red") +
#   geom_vline(xintercept = 1, linetype = "dotted") +
#   theme_bw() +
#   labs(
#     x = expression(R[xy] ~ "(population 1 / population 2)"),
#     y = ""
#   ) +
#   theme(
#     text = element_text(size = 13),
#     plot.margin = margin(5.5, 130, 5.5, 5.5)
#   )+
#   xlim(0.2, 2.5)
# 
# # Combine plot and horizontal title on the right
# poly_mono_sharedSNP_X_title <- ggdraw() +
#   draw_plot(poly_mono_sharedSNP_X, x = 0, width = 0.9) +  # leave space for title on right
#   draw_label("Polyandry\n
#              vs\n
#              Monandry", x = 0.65, y = 0.5,
#              hjust = 0, vjust = 0.5, size = 12, angle = 0, lineheight = 0.55)
# 
# Rxy_sharedSNP_X <- ggarrange(mono_ga1_sharedSNP_X_title,
#                                      poly_ga1_sharedSNP_X_title,
#                                      poly_mono_sharedSNP_X_title,
#                                      ncol=1)
# 
# ggsave(
#   filename = "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/Rxy/Rxy_sharedSNP_X.png",
#   plot = Rxy_sharedSNP_X,
#   device = "png",
#   units="mm",
#   height=120,
#   width=170,
#   dpi=300
# )
