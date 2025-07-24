library(tidyverse)


setwd("~/dispersal/baypass/R/POD")



# arguments supplied with the call to this script
# 1: NPOP
# 2: Omega run XtX output
# 3: SNP position file of XtX run
# 4: POD run XtX output

args <- commandArgs(trailingOnly = TRUE)

npop <- args[1]
omega_xtx_filename <- args[2]
snp_pos_filename <- args[3]
pod_xtx_filename <- args[4]



# load xtx output file from omega run, whichever file supplied to this R script in the sub script
output_file_location <- "~/dispersal/baypass/omega_matrices/"
omega_xtx_file <- paste(output_file_location, omega_xtx_filename, sep="")
xtx <- read_delim(omega_xtx_file,
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
      `xtx_log10(1/pval)` = as.numeric(`log10(1/pval)`)
    )



# load position file from omega run, whichever position file supplied to this R script in the sub script
snp_pos_file_location <- "~/dispersal/baypass/positions_of_snps_used_generate_omega_matrices/"
snp_pos_file <- paste(snp_pos_file_location, snp_pos_filename, sep="")
snp_pos <- read_delim( snp_pos_file, delim = "\t", col_names = F)
colnames(snp_pos) <- c("CHROM", "POS")



# load POD xtx output file using omega run genofile, for whichever position file supplied to this R script in the sub script
pod_xtx_file_location <- "~/dispersal/baypass/simulated_genofiles_omega_runs/POD_run_on_simulated_genofiles_omega_runs/"
pod_xtx_file <- paste(pod_xtx_file_location, pod_xtx_filename, sep="")
pod_xtx <- read.table( pod_xtx_file, header = T) %>%
  mutate( `log10(1/pval)` = log10.1.pval. )
pod_xtx_threshold <- quantile(pod_xtx$M_XtX, 0.995)
pod_xtxst_threshold <- quantile(pod_xtx$XtXst, 0.995)
pod_xtx_p_threshold <- pod_xtx %>%
  filter(XtXst>=pod_xtxst_threshold) %>%
  arrange(desc(XtXst)) %>%
  head(1)%>%
  pull(`log10(1/pval)`)





# add chrom & pos columns from snp position file
# add a cols for whether or not the stats are above the thresholds
xtx <- xtx %>%
  mutate(CHROM = as.factor(snp_pos$CHROM), 
        POS = snp_pos$POS,
        mxtx_greaterthan_POD = as.factor(case_when(   M_xtx >  pod_xtx_threshold ~ 1,
                                            M_xtx <= pod_xtx_threshold ~ 0)),
        xtxst_greaterthan_POD = as.factor(case_when(  xtxst >  pod_xtxst_threshold ~ 1,
                                            xtxst <= pod_xtxst_threshold ~ 0)),
        xtx_p_greaterthan_POD = as.factor(case_when( `xtx_log10(1/pval)` > pod_xtx_p_threshold ~ 1,
                                           `xtx_log10(1/pval)` <= pod_xtx_p_threshold ~ 0))
      )
  
  

  
  
  
# TURN THE NON-CHROM CONTIGS INTO A CONTIGUOUS THING FOR PLOTTING
if (any( grepl("^KQ", xtx$CHROM))){
           
           chrom <- xtx %>% 
            filter(str_detect(CHROM, "^LG"))
           
           non_chrom <- xtx %>% 
             filter(str_detect(CHROM, "^KQ")) %>%
             mutate(POS = cumsum(POS),
                    CHROM = "Non-chrom contigs"
             )
           
           xtx <- bind_rows(non_chrom, chrom)
  }
  
  

########### PLOT XtX manhattan
xtx_plot <- xtx %>%
  ggplot(aes(x = POS, y = M_xtx, colour=mxtx_greaterthan_POD)) +
    geom_point(size = 1.5, alpha = 0.5, stroke = 0) +
    facet_grid(. ~ CHROM, scales = "free_x", switch = "x", space = "free_x") +
    xlab("Chromsome") +
    ylab("XtX") +
    ggtitle(paste("XtX", npop, snp_pos_filename )) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.1, "cm"),
      strip.background = element_blank(),
      strip.placement = "outside",
      legend.position = "none"
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
      scale_colour_manual(values=c("grey60", "red4"))
                   
ggsave(
  filename = paste("ggsave_xtx_", npop , "_", snp_pos_filename, ".png", sep = ""),
  plot = xtx_plot,
  device = "png",
  dpi = 200,
  width = 400,
  height = 120,
  units = "mm"
  )
                   
                   
                   
                   
########### PLOT XtXst
xtxst_plot <- xtx %>%
  ggplot(aes(x = POS, y = xtxst, colour=xtxst_greaterthan_POD)) +
    geom_point( size = 1.5, alpha = 0.5, stroke = 0 ) +
    facet_grid(. ~ CHROM, scales = "free_x", switch = "x", space = "free_x") +
    xlab("Chromsome") +
    ylab("XtXst") +
    ggtitle( paste("XtXst", npop, snp_pos_filename)) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.1, "cm"),
      strip.background = element_blank(),
      strip.placement = "outside",
      legend.position = "none"
      ) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
    scale_colour_manual(values=c("grey60", "red4"))

                  
ggsave(
  filename = paste("ggsave_xtxst_", npop, "_", snp_pos_filename, ".png", sep = ""),
  plot = xtxst_plot,
  device = "png",
  dpi = 300,
  width = 300,
  height = 120,
  units = "mm"
  )
                       
                       
                       
                       
                       
########### PLOT XtX -log10(p) manhattan for the ith run of runs
xtx_p_plot <- xtx %>%
  ggplot(aes(x = POS, y = `xtx_log10(1/pval)`, colour=xtx_p_greaterthan_POD)) +
    geom_point(
      size = 1.5,
      alpha = 0.5,
      stroke = 0
      ) +
    facet_grid( . ~ CHROM, scales = "free_x", switch = "x", space = "free_x") +
    xlab("Chromsome") +
    ylab("XtX -log10(1/p)") +
    ggtitle( paste("XtX p value", npop, snp_pos_filename )) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing = unit(0.1, "cm"),
      strip.background = element_blank(),
      strip.placement = "outside",
      legend.position = "none"
      ) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
      scale_colour_manual(values=c("grey60", "red4"))

                           
ggsave(
  filename = paste("ggsave_xtx_pval_", npop, "_", snp_pos_filename, ".png", sep = ""),
  plot = xtx_p_plot,
  device = "png",
  dpi = 300,
  width = 300,
  height = 120,
  units = "mm"
  )
                           
                  

