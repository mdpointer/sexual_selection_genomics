print(.libPaths())


user_lib <- "~/R/x86_64-pc-linux-gnu-library/4.4"
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE)
}
.libPaths(user_lib)  # Set the library path for this session


# Check if required packages are installed, and install them if not
required_packages <- c("tidyverse", "qvalue")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}


# library(tidyverse)
# library(qvalue)


setwd("~/ss/baypass/R/plots")


# arguments supplied with the call to this script
# 1: NPOP
# 2: Data run XtX output
# 2: Data run BF output
# 2: Data run C2 output
# 5: SNP position file of XtX run
# 6: POD run on sim genofile XtX output

args <- commandArgs(trailingOnly = TRUE)

npop <- as.numeric(args[1])
xtx_filename <- args[2]
BF_filename <- args[3]
C2_filename <- args[4]
snp_pos_filename <- args[5]
pod_xtx_filename <- args[6]
pod_c2_filename <- args[7]

print("objects assigned from args")


# load xtx output file from data run, whichever file supplied to this R script in the sub script
output_file_location <- "~/ss/baypass/data_runs_C2_BF/combined_sex_X/"
xtx_file <- paste(output_file_location, xtx_filename, sep="")
xtx <- read_delim(xtx_file,
                    delim = "\t",
                    col_names = F,
                    skip = 1) %>%
    mutate(X1 = str_squish(X1)) %>%
    separate_wider_delim(
      cols = X1,
      delim = " ",
      names = c(
        "MRK",
        "M_P",
        "SD_P",
        "M_xtx",
        "SD_xtx",
        "xtxst",
        "log10(1/pval)"
      )
    ) %>%
    mutate(
      line_no = 1:nrow(.),
      MRK = as.numeric(MRK),
      M_P = as.numeric(M_P),
      SD_P = as.numeric(SD_P),
      M_xtx = as.numeric(M_xtx),
      SD_xtx = as.numeric(SD_xtx),
      xtxst = as.numeric(xtxst),
      xtxst_recomp_p = 1-2*abs(pchisq(xtxst,df=npop)-0.5),
      xtxst_recomp_neglogp = -log10(xtxst_recomp_p),
      `xtxst_log10(1/pval)` = as.numeric(`log10(1/pval)`),
      xtxst_backcomp_p = 10^(-`xtxst_log10(1/pval)`)
    )
print("xtx file loaded")


# load file to get BF(dB)
output_file_location <- "~/ss/baypass/data_runs_C2_BF/combined_sex_X/"
BF_file <- paste(output_file_location, BF_filename, sep="")
BF <- read_delim(BF_file, delim="\t", col_names = F, skip=1) %>% 
  mutate( X1 = str_squish(X1)) %>% 
  separate_wider_delim(cols=X1, delim=" ", names= c("COVARIABLE", "MRK", "M_Beta", "SD_Beta", "PIP", "BF(dB)")) %>% 
  mutate(line_no = 1:nrow(.),
         MRK = as.numeric(MRK),
         M_Beta = as.numeric(M_Beta),
         SD_Beta = as.numeric(SD_Beta),
         PIP = as.numeric(PIP),
         `BF(dB)` = as.numeric(`BF(dB)`))
print("BF file loaded")



# load file to get C2
output_file_location <- "~/ss/baypass/data_runs_C2_BF/combined_sex_X/"
C2_file <- paste(output_file_location, C2_filename, sep="")
C2 <- read_delim(C2_file, delim="\t", col_names = F, skip=1) %>%
  mutate( X1 = str_squish(X1)) %>%
  separate_wider_delim(cols=X1, delim=" ", names= c("CONTRAST", "MRK", "M_C2", "SD_C2", "C2_std", "log10(1/pval)")) %>%
  mutate(line_no = 1:nrow(.),
         CONTRAST = as.numeric(CONTRAST),
         MRK = as.numeric(MRK),
         M_C2 = as.numeric(M_C2),
         SD_C2 = as.numeric(SD_C2),
         C2_std = as.numeric(C2_std),
         C2_recomp_p = pchisq(C2_std, df=1,lower.tail=F),
         C2_recomp_neglogp = -log10(C2_recomp_p),
         `C2_log10(1/pval)` = as.numeric(`log10(1/pval)`),
         C2_backcomp_p = 10^(-`C2_log10(1/pval)`))
print("C2 file loaded")


# load position file from omega run, whichever position file supplied to this R script in the sub script
snp_pos_file_location <- "~/ss/baypass/positions_of_snps_used_baypass_data_runs/"
snp_pos_file <- paste(snp_pos_file_location, snp_pos_filename, sep="")
snp_pos <- read_delim( snp_pos_file, delim = "\t", col_names = F)
colnames(snp_pos) <- c("CHROM", "POS")
print("snps file loaded")

# 
# # load POD xtx output file using omega run genofile, for whichever position file supplied to this R script in the sub script
quantile_to_use_for_POD <- 0.999
pod_xtx_file_location <- "~/ss/baypass/POD/simulated_genofiles_data_runs/combined_sex_X/"
pod_xtx_file <- paste(pod_xtx_file_location, pod_xtx_filename, sep="")
pod_xtx <- read.table( pod_xtx_file, header = T) %>%
  mutate( `log10(1/pval)` = log10.1.pval. )
pod_xtx_threshold <- quantile(pod_xtx$M_XtX, quantile_to_use_for_POD)
pod_xtxst_threshold <- quantile(pod_xtx$XtXst, quantile_to_use_for_POD)
pod_xtx_p_threshold <- pod_xtx %>%
  filter(XtXst>=pod_xtxst_threshold) %>%
  arrange(desc(XtXst)) %>%
  head(1)%>%
  pull(`log10(1/pval)`)
cat("The ", quantile_to_use_for_POD, " xtx quantile is: ", pod_xtx_threshold, "\n")


# 
# # load POD C2 output file
pod_c2_file_location <- "~/ss/baypass/POD/simulated_genofiles_data_runs/combined_sex_X/"
pod_c2_file <- paste(pod_c2_file_location, pod_c2_filename, sep="")
pod_c2 <- read.table( pod_c2_file, header = T)
pod_c2_threshold_0.999 <- quantile(pod_c2$C2_std, 0.999)

cat("The 0.999 c2 quantile is: ", pod_c2_threshold_0.999, "\n")




# add chrom & pos columns from snp position file
# add a cols for whether or not the stats are above the thresholds
xtx <- xtx %>%
  mutate(CHROM = as.factor(snp_pos$CHROM), 
        POS = snp_pos$POS,
        mxtx_greaterthan_POD = as.factor(case_when(   M_xtx >  pod_xtx_threshold ~ 1,
                                            M_xtx <= pod_xtx_threshold ~ 0)),
        xtxst_greaterthan_POD = as.factor(case_when(  xtxst >  pod_xtxst_threshold ~ 1,
                                            xtxst <= pod_xtxst_threshold ~ 0)),
        xtx_p_greaterthan_POD = as.factor(case_when( `xtxst_log10(1/pval)` > pod_xtx_p_threshold ~ 1,
                                           `xtxst_log10(1/pval)` <= pod_xtx_p_threshold ~ 0))
      ) %>% 
  mutate(CHROM_POS = paste(CHROM, POS, sep="_"))

BF <- BF %>%
  mutate(CHROM = as.factor(snp_pos$CHROM), 
         POS = snp_pos$POS) %>% 
  mutate(CHROM_POS = paste(CHROM, POS, sep="_"))

C2 <- C2 %>%
  mutate(CHROM = as.factor(snp_pos$CHROM), 
         POS = snp_pos$POS ) %>% 
  mutate(CHROM_POS = paste(CHROM, POS, sep="_"))


# join the dataframes together by the CHROM_POS col
xtx_BF <- full_join(xtx, BF, by="CHROM_POS")
dd <- full_join(xtx_BF, C2, by="CHROM_POS")                     ## from here the full dataset is called dd
# remove intermediate dataframes
rm(list= "xtx", "BF", "C2", "xtx_BF")
print("data processed and joined by snp")




# ###
# ####  Mxtx window analysis
# ###
# 
# window_size <- 100000
# window_size_kb <- window_size/1000
# step_size = window_size/2 # DO NOT CHANGE THIS LINE, THE COMPUTATION IS BASED ON THE STEP BEING HALF THE WINDOW WIDTH
# n_window_threshold = 2
# 
# # calc counts in sliding windows
# window_mxtx_counts <- dd %>% 
#   mutate(window1 = (POS %/% window_size) + 1,
#          window2 = ((POS - window_size/2) %/% window_size)+1.5) %>%
#   pivot_longer(cols = c(window2, window1),
#                names_to = "window_type",
#                values_to = "window"
#   ) %>%
#   group_by(CHROM, window) %>% 
#   summarise( count = sum(M_xtx > pod_xtx_threshold)) %>% 
#   mutate(startPOS = (window * window_size) - window_size,
#          endPOS   = (window * window_size) -1) %>% 
#   filter(count >= n_window_threshold)
# 
# windows_with_gteq_2_snps_above_POD_threshold <- window_mxtx_counts %>% nrow()
# cat("There are ", windows_with_gteq_2_snps_above_POD_threshold, window_size_kb,"kb windows with >=2 snps in excess of the ", quantile_to_use_for_POD, " M_XtX threshold", "\n")



###
####  C2 significance analysis - SNPs
###

pod_c2_threshold <- pod_c2_threshold_0.999 ## insert the POD threshold object from above

dd %>%
  filter(C2_std > pod_c2_threshold_0.999) %>%
  write.table(paste("~/ss/baypass/outlier_output/", C2_filename, "baypass_c2_outlier_snps_pod0.999", sep=""),
              quote=F, sep="\t", row.names=F) ## insert quantile in outut filename

dd %>%
  filter(C2_std > pod_c2_threshold_0.999) %>%
  filter(`BF(dB)` >= 20) %>% 
  write.table(paste("~/ss/baypass/outlier_output/", C2_filename, "baypass_c2BF_outlier_snps_pod0.999_BF20", sep=""),
              quote=F, sep="\t", row.names=F) ## insert quantile in outut filename





####
### PLOTS

# ### validation plots
# 
# # check that my recomputed XTXst P is the same as back-computing p from log10(1/p)
# xtxst_recompP_vs_backcompP <- dd %>%
#   ggplot(aes(x=xtxst_recomp_p, y=xtxst_backcomp_p)) +
#   geom_point(alpha=0.5) +
#   theme_bw() +
#   ggtitle( paste("Check_computation_of_xtxst_pvalues", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_"))
# 
# ggsave(
#   filename = paste("ggsave_xtxst_check_pval_comp_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
#   plot = xtxst_recompP_vs_backcompP,
#   device = "png",
#   dpi = 300,
#   width = 200,
#   height = 120,
#   units = "mm"
# )
# 
# 
# # check that my recomputed C2 P is the same as back-computing p from log10(1/p)
# C2_recompP_vs_backcompP <- dd %>%
#   ggplot(aes(x=C2_recomp_p, y=C2_backcomp_p)) +
#   geom_point(alpha=0.5) +
#   theme_bw() +
#   ggtitle( paste("Check_computation_of_C2_pvalues", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_"))
# 
# ggsave(
#   filename = paste("ggsave_C2_check_pval_comp_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
#   plot = C2_recompP_vs_backcompP ,
#   device = "png",
#   dpi = 300,
#   width = 200,
#   height = 120,
#   units = "mm"
# )
# print("xtxst and C2 Pval checks plotted")
# 
# 
# 
# ########### PLOT XtX histogram
# xtx_hist_plot <- dd %>%
#   ggplot(aes(x=M_xtx)) +
#   geom_histogram(bins=20, boundary=0, closed="left") +
#   theme_bw() +
#   ggtitle( paste("xtx_hist", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_"))
# 
# ggsave(
#   filename = paste("ggsave_xtx_hist_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
#   plot = xtx_hist_plot,
#   device = "png",
#   dpi = 300,
#   width = 200,
#   height = 120,
#   units = "mm"
# )
# print("xtx hist plotted")
# 
# 
# ########### PLOT XtXst histogram
# xtxst_hist_plot <- dd %>%
#   ggplot(aes(x=xtxst)) +
#   geom_histogram(bins=20, boundary=0, closed="left") +
#   theme_bw() +
#   ggtitle( paste("xtxst_hist", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_"))
# 
# ggsave(
#   filename = paste("ggsave_xtxst_hist_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
#   plot = xtxst_hist_plot,
#   device = "png",
#   dpi = 300,
#   width = 200,
#   height = 120,
#   units = "mm"
# )
# print("xtxst hist plotted")
# 
# ########### PLOT XtXst pval histogram
# xtxst_p_hist_plot <- dd %>%
#   ggplot(aes(x=xtxst_recomp_p)) +
#   geom_histogram(binwidth =0.05, boundary=0, closed="left") +
#   theme_bw() +
#   ggtitle( paste("xtxst_pval_hist", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_"))
# 
# ggsave(
#   filename = paste("ggsave_xtxst_pval_hist_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = "" ),
#   plot = xtxst_p_hist_plot,
#   device = "png",
#   dpi = 300,
#   width = 200,
#   height = 120,
#   units = "mm"
# )
# print("xtxst pval hist plotted")
# 
# ########### PLOT C2_std histogram
# C2_std_plot <- dd %>%
#   ggplot(aes(x=C2_std)) +
#   geom_histogram(binwidth =0.05, boundary=0, closed="left") +
#   theme_bw() +
#   ggtitle( paste("C2_hist", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_"))
# 
# ggsave(
#   filename = paste("ggsave_C2_std_hist_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
#   plot = C2_std_plot,
#   device = "png",
#   dpi = 300,
#   width = 200,
#   height = 120,
#   units = "mm"
# )
# print("C2std hist plotted")
# 
# ########### PLOT C2 pval histogram
# C2_p_hist_plot <- dd %>%
#   ggplot(aes(x=C2_recomp_p)) +
#   geom_histogram(binwidth =0.05, boundary=0, closed="left") +
#   theme_bw() +
#   ggtitle( paste("C2_pval_hist", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_"))
# 
# ggsave(
#   filename = paste("ggsave_C2_pval_hist_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
#   plot = C2_p_hist_plot,
#   device = "png",
#   dpi = 300,
#   width = 200,
#   height = 120,
#   units = "mm"
# )
# print("C2 pval hist plotted")




### XTX

# ########### PLOT XtX manhattan
# Mxtx_plot <-
#   ggplot() +
#   geom_rect(data=window_mxtx_counts,
#             aes(xmin=startPOS,xmax=endPOS,ymin=-Inf,ymax=Inf),
#             alpha=0.3, fill="blue4") +
#     geom_point(data=dd,
#                aes(x = POS, y = M_xtx, colour=mxtx_greaterthan_POD),
#                size = 1.5, alpha = 0.5, stroke = 0) +
#     geom_hline(yintercept = pod_xtx_threshold, color = "red4", linetype = "dashed") +
#     facet_grid(. ~ CHROM, scales = "free_x", switch = "x", space = "free_x") +
#     xlab("Chromsome") +
#     ylab("M_XtX") +
#     ggtitle(paste("M_XtX", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_")) +
#     theme_classic() +
#     theme(
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       panel.spacing = unit(0.1, "cm"),
#       strip.background = element_blank(),
#       strip.placement = "outside",
#       legend.position = "none"
#       ) +
#       scale_x_continuous(expand = c(0, 0)) +
#       # scale_y_continuous( limits = c(0, NA))expand = c(0, 0), #+
#       scale_colour_manual(values=c("grey60", "red4"))
# 
# ggsave(
#   filename = paste("ggsave_Mxtx_", substr(xtx_filename, 19, nchar(xtx_filename)-24), "_", quantile_to_use_for_POD, "quant", "_", window_size_kb, "windows", ".png", sep = ""),
#   plot = Mxtx_plot,
#   device = "png",
#   dpi = 200,
#   width = 400,
#   height = 120,
#   units = "mm"
#   )
# print("xtx manhattan plotted")




# ########### PLOT XtXst manhattan
# xtxst_plot <- dd %>%
#   ggplot(aes(x = POS, y = xtxst)) +
#     geom_point( size = 1.5, alpha = 0.5, stroke = 0, colour="grey60") +
#     facet_grid(. ~ CHROM, scales = "free_x", switch = "x", space = "free_x") +
#     xlab("Chromsome") +
#     ylab("XtXst") +
#     ggtitle( paste("XtXst", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_")) +
#     theme_classic() +
#     theme(
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       panel.spacing = unit(0.1, "cm"),
#       strip.background = element_blank(),
#       strip.placement = "outside",
#       legend.position = "none"
#       ) +
#     scale_x_continuous(expand = c(0, 0))
#     # scale_y_continuous(limits = c(0, NA)) expand = c(0, 0),  #+
#     #scale_colour_manual(values=c("grey60")) #, "red4"))
# 
# ggsave(
#   filename = paste("ggsave_xtxst_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
#   plot = xtxst_plot,
#   device = "png",
#   dpi = 300,
#   width = 300,
#   height = 120,
#   units = "mm"
#   )
# print("xtxst manhattan plotted")




# 
# ### BF
# 
# ########### PLOT BF manhattan
# BF_plot <- dd %>%
#   ggplot(aes(x = POS, y = `BF(dB)`)) +
#   geom_point(size = 1.5, alpha = 0.5, stroke = 0, colour="grey60") +
#   facet_grid( . ~ CHROM, scales = "free_x", switch = "x", space = "free_x") +
#   xlab("Chromsome") +
#   ylab("BF") +
#   ggtitle( paste("BF", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_")) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     panel.spacing = unit(0.1, "cm"),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     legend.position = "none"
#   ) +
#   scale_x_continuous(expand = c(0, 0))
# #   scale_y_continuous(limits = c(NA, NA) ), expand = c(0, 0)) #+
#   #scale_colour_manual(values=c("grey60")) #, "red4"))
# 
# ggsave(
#   filename = paste("BF_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
#   plot = BF_plot,
#   device = "png",
#   dpi = 300,
#   width = 300,
#   height = 120,
#   units = "mm"
# )
# print("BF manhattan plotted")




########### PLOT C2 manhattan
C2_plot <- dd %>% 
  filter(str_detect(CHROM, "^LG")) %>% 
  mutate(CHROM = fct_relevel(CHROM, "LGX"),
         point_col = 1,
         point_col = case_when(C2_std > pod_c2_threshold ~ 2, # assign a highlight colour to points above the sig threshold
                               C2_std <= pod_c2_threshold ~ point_col)) %>%
  ggplot(aes(x = POS, y = C2_std, col= as.factor(point_col), size=as.factor(point_col))) +
  geom_point(alpha = 0.5, stroke = 0) +
  facet_grid( . ~ CHROM, scales = "free_x", switch = "x", space = "free_x") +
  xlab("Genome position") +
  ylab("C2") +
  # ggtitle( paste("C2", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_")) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing = unit(0.1, "cm"),
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.position = "none",
    # plot.title = element_text(face = "bold"),
    text = element_text(size = 20)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(
    breaks = c(0.0, 2.5, 5.0, 7.5),
    limits = c(0, NA),
    expand = c(0, 0)
  )+
  scale_colour_manual(values=c("grey60", "red4")) +
  scale_size_manual(values=c(1.5, 1.5, 2.5)) +
    geom_hline(yintercept = pod_c2_threshold_0.999, color = "red4", linetype = "dashed")

ggsave(
  filename = paste("ggsave_C2_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
  plot = C2_plot,
  device = "png",
  dpi = 300,
  width = 55,
  height = 120,
  units = "mm"
)
cat("C2 manhattan plotted \n")

ggsave(
  filename = paste("ggsave_C2_", substr(xtx_filename, 19, nchar(xtx_filename)-24), "_wide.png", sep = ""),
  plot = C2_plot,
  device = "png",
  dpi = 300,
  width = 400,
  height = 120,
  units = "mm"
)
cat("C2 manhattan plotted \n")







# ####
# ####
# #### Q VALUES
# 
# 
# #### COMPUTE Q VALUES
# # use qvalue package to convert pvalues to qvalues & add to dataframe
# xtxst_recomp_p <- dd$xtxst_recomp_p
# C2_recomp_p <- dd$C2_recomp_p
# xtxst_q <- qvalue(p=xtxst_recomp_p)
# C2_q <- qvalue(p=C2_recomp_p)
# rm(list= "xtxst_recomp_p", "C2_recomp_p")
# 
# dd <- dd %>%
#   mutate( xtxst_qval = xtxst_q$qvalues,
#           C2_qval = C2_q$qvalues,
#           `xtxst_log10(1/qval)` = -log10(xtxst_q$qvalues),
#           `C2_log10(1/qval)` = -log10(C2_q$qvalues)
#   )
# print("Qvalues computed")
# 
# 
# 
# 
# ########### PLOT XtX qval histogram
# xtxst_q_hist_plot <- dd %>%
#   ggplot(aes(x=xtxst_qval)) +
#   geom_histogram(binwidth =0.05, boundary=0, closed="left") +
#   theme_bw() +
#   ggtitle(paste("xtx_qval_hist", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_"))
# 
# 
# ggsave(
#   filename = paste("ggsave_xtxst_qval_hist_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
#   plot = xtxst_q_hist_plot,
#   device = "png",
#   dpi = 300,
#   width = 200,
#   height = 120,
#   units = "mm"
# )
# print("xtxst qval hist plotted")
# 
# 
# ########### PLOT C2 qval histogram
# C2_q_hist_plot <- dd %>%
#   ggplot(aes(x=C2_qval)) +
#   geom_histogram(binwidth =0.05, boundary=0, closed="left") +
#   theme_bw() +
#   ggtitle( paste("C2_qval_hist", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_"))
# 
# ggsave(
#   filename = paste("ggsave_C2_qval_hist_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
#   plot = C2_q_hist_plot,
#   device = "png",
#   dpi = 300,
#   width = 200,
#   height = 120,
#   units = "mm"
# )
# print("C2 qval hist plotted")




########### PLOT BF vs C2

plot_bf_C2 <- dd  %>%
  mutate(bf_threshold = case_when(`BF(dB)` > 20 ~ 1,
                                  TRUE ~ 0),
         C2_threshold = case_when(C2_std > pod_c2_threshold_0.999 ~ 1,
                                  TRUE ~ 0),
         bfC2_threshold = case_when(bf_threshold ==1 & C2_threshold == 1 ~ 1,
                                  TRUE ~ 0)) %>% 
  ggplot(aes(x=`BF(dB)`, y=C2_std, col=as.factor(bfC2_threshold))) +
  geom_hline(yintercept = pod_c2_threshold_0.999, color = "red4", linetype = "dashed") +
  geom_vline(xintercept = 20, colour="red4", linetype="dashed") +
  geom_point(alpha=0.5) +
  theme_bw() +
  scale_colour_manual(values=c("black", "red4"))
  ggtitle( paste("BF_vs_C2", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_"))


ggsave(filename= paste("BF_v_C2", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep=""),
       plot= plot_bf_C2,
       device= "png",
       dpi= 300,
       width=250,
       height=250,
       units= "mm")
print("BF vs C2 plotted")



# ########### PLOT XtX -log10(1/q) manhattan
# xtxst_q_plot <- dd %>%
#   ggplot(aes(x = POS, y = `xtxst_log10(1/qval)`)) +
#   geom_point(size = 1.5, alpha = 0.5, stroke = 0, colour="grey60") +
#   facet_grid( . ~ CHROM, scales = "free_x", switch = "x", space = "free_x") +
#   xlab("Chromsome") +
#   ylab("XtXst -log10(qvalue)") +
#   ggtitle( paste("XtXst_q-value", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_" )) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     panel.spacing = unit(0.1, "cm"),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     legend.position = "none"
#   ) +
#   scale_x_continuous(expand = c(0, 0))
# #   scale_y_continuous(limits = c(0, NA)) expand = c(0, 0), ) #+
#   #scale_colour_manual(values=c("grey60")) #, "red4"))
# 
# ggsave(
#   filename = paste("ggsave_xtxst_qval_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
#   plot = xtxst_q_plot,
#   device = "png",
#   dpi = 300,
#   width = 300,
#   height = 120,
#   units = "mm"
# )
# print("xtxst Qval manhattan plotted")
# 
# 
# 
# ### C2
# 
# ########### PLOT C2 -log10(1/q) manhattan
# C2_q_plot <- dd %>%
#   ggplot(aes(x = POS, y = `C2_log10(1/qval)`)) +
#   geom_point( size = 1.5, alpha = 0.5, stroke = 0, colour="grey60") +
#   facet_grid( . ~ CHROM, scales = "free_x", switch = "x", space = "free_x") +
#   xlab("Chromsome") +
#   ylab("C2 -log10(qvalue)") +
#   ggtitle( paste("C2_q-value", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_")) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_blank(),
#     axis.ticks.x = element_blank(),
#     panel.spacing = unit(0.1, "cm"),
#     strip.background = element_blank(),
#     strip.placement = "outside",
#     legend.position = "none"
#   ) +
#   scale_x_continuous(expand = c(0, 0))
#   # scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
#   #scale_colour_manual(values=c("grey60")) #, "red4"))
# 
# ggsave(
#   filename = paste("ggsave_C2_qval_", substr(xtx_filename, 19, nchar(xtx_filename)-24), ".png", sep = ""),
#   plot = C2_q_plot,
#   device = "png",
#   dpi = 300,
#   width = 300,
#   height = 120,
#   units = "mm"
# )
# cat("C2 Qval manhattan plotted \n")






# ##### PLOT C2 ZOOMED IN FOR A A SELECTION OF INDIVIDUAL PEAKS
# 
# ########### PLOT C2 manhattan
# local_C2_plot <- dd %>% 
#   filter(str_detect(CHROM, "^LG")) %>% 
#   filter(CHROM=="LG2" & POS>=9489243 & POS<=10499733)
#   mutate(point_col = case_when(C2_std > pod_c2_threshold ~ 2, # assign a highlight colour to points above the sig threshold
#                                C2_std <= pod_c2_threshold ~ 1)) %>%
#   ggplot(aes(x = POS, y = C2_std, col= as.factor(point_col))) +
#   geom_point(alpha = 0.5, stroke = 0, size=2) +
#     geom_line() +
#   xlab("Genome position") +
#   ylab("C2") +
#   # ggtitle( paste("C2", substr(xtx_filename, 19, nchar(xtx_filename)-24), sep="_")) +
#   theme_classic() +
#   theme(
#     axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#     legend.position = "none",
#     # plot.title = element_text(face = "bold"),
#     text = element_text(size = 20)
#   ) +
#   scale_x_continuous(expand = c(0, 0)) +
#   # scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
#   scale_colour_manual(values=c("grey30", "red4")) +
#   geom_hline(yintercept = pod_c2_threshold_0.999, color = "red4", linetype = "dashed")
# 
# ggsave(
#   filename = "LG2_FST_around_c2SNP_loc9994244.png",
#   plot = local_C2_plot,
#   device = "png",
#   dpi = 300,
#   width = 100,
#   height = 80,
#   units = "mm"
# )
# cat("C2 manhattan plotted \n")


cat("End of script")



