library(tidyverse)
library(ggpubr)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggh4x)


## read in population file to add treatment info
# there is one 'wrong' name in this list, but it is also wrong in the VCF so leave it here and correct after
source("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/sample_info/generate_population.list.full.R")
population_info <- generate_ss_population.info.full.list("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/scripts/sample_info/sample_names_file_182_from_VCF.txt") %>%
  rename(ind = name) %>%
  mutate( line_rep = case_when( ind== "polym4" ~ "C",
                                TRUE ~ line_rep)) %>%
  mutate(treat_rep = paste(treatment, line_rep, sep="_"))



##########
##########
#########

data_file <- "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/plink_ROH/plink_ROH_run3.hom.indiv"

dd_roh <- read.delim(data_file, header=TRUE, sep = "") %>%
  mutate(ind = case_when(IID == "L004" ~ paste(FID, IID, sep="_"),
                         TRUE ~ FID)) %>% 
  mutate(treatment = substr(ind, 1, 2)) %>% 
  mutate(treatment = case_when( treatment == "GA" ~ "GA1",
                                treatment == "Mb" ~ "male biased",
                                treatment == "Fb" ~ "female biased",
                                treatment == "mo" ~ "mono",
                                treatment == "po" ~ "poly")) %>% 
  # filter(treatment %in% c("GA1","mono","poly")) %>% 
  select(-FID, -IID, -PHE) %>% 
  mutate(Froh = KB / (147631796 / 1000))


dd <- full_join(dd_roh, population_info) %>% 
  select(-seq_run) %>% 
  filter(treatment %in% c("mono","poly"))

Froh_plot <- dd %>% 
  mutate(treatment = case_when( treatment == "mono" ~ "Monandry",
                                treatment == "poly" ~ "Polyandry")) %>% 
  rename(`SS regime` = treatment,
         `Replicate line` = line_rep) %>% 
  ggplot( aes( x= `SS regime`, y=Froh, fill=`SS regime`)) +
  geom_jitter(aes(shape=`Replicate line`),size=4,position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), alpha = 0.4) +  
  geom_boxplot(alpha=0.5, outlier.shape=NA)+
  theme_bw() +
  theme(text = element_text(size = 20)) +
  scale_fill_manual(values=c("dodgerblue4","firebrick")) +
  labs(
    # title="Nonsense",
    y = expression(F[ROH]),
    x="Sexual selection regime") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/ss_Froh.png",
       Froh_plot,
       device="png",
       units="mm",
       width=200,
       height=150)







#################
#################


data_file <- "/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/plink_ROH/plink_ROH_run3.hom"

dd_roh2 <- read.delim(data_file, header=TRUE, sep = "") %>%
  mutate(ind = case_when(IID == "L004" ~ paste(FID, IID, sep="_"),
                         TRUE ~ FID)) %>% 
  mutate(treatment = substr(ind, 1, 2)) %>% 
  mutate(treatment = case_when( treatment == "GA" ~ "GA1",
                                treatment == "Mb" ~ "male biased",
                                treatment == "Fb" ~ "female biased",
                                treatment == "mo" ~ "mono",
                                treatment == "po" ~ "poly")) %>% 
  # filter(treatment %in% c("GA1","mono","poly")) %>% 
  select(-FID, -IID, -PHE) %>% 
  mutate(Froh = KB / (147631796 / 1000))



dd <- full_join(dd_roh2, population_info) %>% 
  select(-seq_run) %>% 
  filter(treatment %in% c("mono","poly"))




# Define genome length in kb (replace with your actual genome length)
total_genome_length_kb <- (147631796 / 1000)

# Create length_class factor as before
dd$length_class <- cut(dd$KB / 1000,
                       breaks = c(0, 0.4, 0.6, 0.8, 1, Inf),
                       labels = c("<0.4 Mb", "0.4–0.6 Mb", "0.6–0.8 Mb", "0.8-1 Mb", ">1"))

# Summarise total ROH length (in kb) per individual, treatment and length class
summary <- dd %>%
  group_by(treatment, ind, length_class) %>%
  summarise(total_length_kb = sum(KB), .groups = "drop") %>%
  mutate(prop_genome = total_length_kb / total_genome_length_kb)

# Calculate total ROH length per individual to order individuals within treatments
total_roh <- summary %>%
  group_by(treatment, ind) %>%
  summarise(total_length = sum(total_length_kb), .groups = "drop") %>%
  arrange(treatment, desc(total_length))

# Reorder 'ind' factor by descending total ROH length within treatment
summary$ind <- factor(summary$ind, levels = total_roh$ind)

# Treatment colors for facet strip backgrounds
treatment_colors <- c(
  "TreatmentA" = alpha("dodgerblue4", 0.5),
  "TreatmentB" = alpha("firebrick", 0.5)
)

# Bar fill colors for length classes (fix empty colors)
length_class_colors <- c(
  "<0.4 Mb" = "gold",
  "0.4–0.6 Mb" = "darkorange3",
  "0.6–0.8 Mb" = "red4",
  "0.8-1 Mb" = "black",
  ">1" = "black"
)

custom_labels <- c(
  "mono" = "Monandry",
  "poly" = "Polyandry"
)

# Plot
Froh_dist_plot <- ggplot(summary, aes(x = ind, y = prop_genome, fill = length_class)) +
  geom_bar(stat = "identity") +
  labs(x = "Individual", y = "Prop of genome in ROHs", fill = "ROH length class") +
  scale_fill_manual(values = length_class_colors) +
  theme_bw() +
  facet_wrap2(~treatment, scales = "free_x", nrow = 1,
              labeller = labeller(treatment = custom_labels),
              strip = strip_themed(
                background_x = elem_list_rect(
                  fill = treatment_colors,
                  color = "black"
                ),
                text_x = elem_list_text(
                )
              )) +
  theme(text = element_text(size = 20),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom"
  )

ggsave("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/ss_Froh_dist.png",
       Froh_dist_plot,
       device="png",
       units="mm",
       width=200,
       height=150)

combi <- ggarrange( Froh_plot,
                    Froh_dist_plot,
                    nrow=2,
                    labels="AUTO",
                    font.label = list(size = 26, face = "bold", color = "black"))


ggsave("/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/ss_Froh_combi.png",
       combi,
       device="png",
       units="mm",
       width=300,
       height=300)

#############
############ model Froh

dd_mod <- dd %>% 
  filter(treatment %in% c("mono","poly"))

ggplot(dd_mod, aes(x = Froh)) +
  geom_histogram(binwidth = 0.001, fill = "steelblue", color = "black") +
  labs(title = "Histogram of Froh", x = "Froh", y = "Count") +
  theme_minimal()



model <- lmer(Froh ~ treatment + (1|treat_rep), data = dd_mod)
summary(model)

qqnorm(resid(model)); qqline(resid(model))
plot(model)  # plots residuals and random effects


dd_mod %>% 
  group_by(treatment) %>% 
  summarise(mean_Froh = mean(Froh),
            se_Froh = sd(Froh) / sqrt(n()))






