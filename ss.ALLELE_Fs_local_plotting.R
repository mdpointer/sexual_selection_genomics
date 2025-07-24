library(tidyverse)

## Local plotting

dd_prop <- read_rds(file="/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/data/allele_Fs/processed_minor_allele_F_proportions_for_AFS_plot.rds")


out_plot <- dd_prop %>%
  #filter(pop %in% c("GA1","mbiasedall","fbiasedall","monoall","polyall")) %>%
  filter(pop %in% c("GA1","mbiasedall","mbiasedA","mbiasedB","mbiasedC","fbiasedall","fbiasedA","fbiasedB","fbiasedC","monoall","monoA","monoB","monoC","polyall","polyA","polyB","polyC")) %>%
  filter(chrom == "autosomes") %>%
  ggplot( aes(x=bin, y=prop, fill=var_cat))+
  geom_bar(stat = "identity", position = "dodge")+
  facet_wrap(~pop, nrow = )+
  theme_bw()+
  scale_fill_manual(values=c("chartreuse4","sienna2","royalblue", "grey50"))+
  labs(
    x="Allele frequency",
    y="Fraction of sites",
    fill="Variant category"
  )+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title.x=element_text(size=20),
    axis.title.y=element_text(size=20),
    strip.text = element_text(size = 14),
    legend.position = "top"  # Place legend at the top
  )

ggsave(plot=out_plot,
       file="/Users/mdp/Google Drive/WORK/tribolium_postdoc/ss_lines/ss_genomics/plots/allele_Fs_per_replicate.png",
       device="png",
       units="mm",
       height=300,
       width=300)




dd_prop %>% 
  filter( !(grepl("all", pop))) %>% 
  mutate(expt = case_when( grepl( "biased", pop) ~ "MbFb",
                           grepl( "mono", pop)   ~ "MonoPoly",
                           grepl( "poly", pop)   ~ "MonoPoly",
                           grepl( "GA1", pop)    ~ "GA1"),
         treatment = case_when( grepl( "fbiased", pop) ~ "female_biased",
                                grepl( "mbiased", pop) ~ "male_biased",
                                grepl( "mono", pop)   ~ "mono",
                                grepl( "poly", pop)   ~ "poly",
                                grepl( "GA1", pop)    ~ "GA1"),
         sex = case_when( grepl("males", pop) ~ "single_sex",
                          TRUE ~ "all_sexes")) %>% 
  filter(chrom == "autosomes",
         sex=="all_sexes") %>% 
  filter(pop=="mbiasedB") %>% 
  filter(var_cat!="missense_tol") %>% 
  ggplot( aes(x=bin, y=prop, group=var_cat, col = var_cat))+
  geom_line(linewidth=1.5)+
  facet_wrap(~pop, nrow = )+
  theme_bw()+
  scale_fill_manual(values=c("chartreuse4","sienna2","royalblue", "grey50"))+
  labs(
    x="Allele frequency",
    y="Proportion of sites",
    fill="Variant category"
  ) +
  theme(text = element_text(size = 20))
  
  
  
dd_prop %>% 
  filter( !(grepl("all", pop))) %>% 
  mutate(expt = case_when( grepl( "biased", pop) ~ "MbFb",
                           grepl( "mono", pop)   ~ "MonoPoly",
                           grepl( "poly", pop)   ~ "MonoPoly",
                           grepl( "GA1", pop)    ~ "GA1"),
         treatment = case_when( grepl( "fbiased", pop) ~ "female_biased",
                                grepl( "mbiased", pop) ~ "male_biased",
                                grepl( "mono", pop)   ~ "mono",
                                grepl( "poly", pop)   ~ "poly",
                                grepl( "GA1", pop)    ~ "GA1"),
         sex = case_when( grepl("males", pop) ~ "single_sex",
                          TRUE ~ "all_sexes")) %>% 
  filter(chrom == "autosomes",
         sex=="all_sexes") %>% 
  group_by(var_cat, bin, expt, treatment) %>% 
  summarise( mean_prop = mean(prop),
             sd_prop = sd(prop),
             n=n()) %>% 
  mutate( se = sd_prop / sqrt(n)) %>% 
  filter(treatment!="GA1",
         var_cat !="missense_tol") %>% 
  filter(treatment=="female_biased") %>% 
  ggplot( aes(x=bin, y=mean_prop, group=var_cat, col = var_cat))+
  geom_line(linewidth=1.5)+
  geom_errorbar(aes(ymin=mean_prop-sd_prop, ymax=mean_prop+sd_prop), width=.2,
                position=position_dodge(0.05)) +
  facet_wrap(~treatment, nrow = )+
  theme_bw()+
  scale_fill_manual(values=c("chartreuse4","sienna2","royalblue", "grey50"))+
  labs(
    x="Allele frequency",
    y="Fraction of sites",
    fill="Variant category"
  )







