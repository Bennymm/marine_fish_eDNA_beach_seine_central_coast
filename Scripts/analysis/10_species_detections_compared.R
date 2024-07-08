#Format data and compare read abundances, biomass densities, and abundances
#create a long data format containing both methods for subsequent analysis

# Packages and data ####
library(tidyverse)
library(lubridate)
library(lme4)
library(scales)
library(dplyr)
library(ggpubr)
library(ggtext)
library(ggside)

#library(plyr) #required for rbind.fill

bs_long1 <- read_rds("Data/2022_10_31/derived_data/bs_long.rds")
eDNA_long1 <- read_rds("Data/2022_10_31/derived_data/eDNA_long.rds")
bs_taxonomy <- read_rds("Data/2022_10_31/derived_data/BS_taxonomy_reconciled.rds")
eDNA_taxonomy <- read_rds("Data/2022_10_31/derived_data/ASV_taxonomy_reconciled.rds")
taxonomy_shared <- read_rds("Data/2022_10_31/derived_data/taxonomy_combined.rds")

# collapse taxonomy and derive read_index ####
#because we are directly comparing metrics of abundance or density for each taxa collapse taxonomy
#to LCT_shared in the taxonomy tables
bs_long <- bs_long1 %>%
  merge(., bs_taxonomy[c("query", "LCT_shared")], by = "query") %>%
  dplyr::select(-c("query")) %>%
  group_by(dat_site, site, date, LCT_shared) %>%
  summarize_at(c("abundance", "weight", "count_per_m2", "count_per_m3", "grams_per_m2", "grams_per_m3"), sum) %>%
  mutate(p_a_bs = if_else(weight == 0, 0, 1)) %>% #recalculate presence / absence after summary
  #rename(p_a_bs = p_a)%>% 
  ungroup()  %>%
  dplyr::select(-c("date", "site"))

e1 <- eDNA_long1 %>%
  merge(., eDNA_taxonomy[c("ASV", "LCT_shared")], by = "ASV") %>%            
  dplyr::select(-c("ASV")) %>%
  group_by(dat_site, site, date, hakai_link_date, LCT_shared) %>%
  summarize_at(c("count"), sum)%>% 
  mutate(p_a_eDNA = if_else(count == 0, 0, 1)) %>% #recalculate presence / absence after summary
  ungroup() 
e2 <- e1 %>%
  dplyr::select(c("dat_site", "LCT_shared", "count"))%>%
  pivot_wider(names_from = LCT_shared, values_from = count) %>%
  column_to_rownames("dat_site")
#calculate read index
PROP <- function(df) {sweep(df, MARGIN=2, colSums(df), FUN="/")}
#Take indices using a table with proportions of reads with families/OTUs in rows, samples in columns
INDEX <- function(df) {sweep(df, MARGIN=1, STATS=apply(df, MARGIN = 1, max), FUN="/")}
eDNA_long <- e2 %>%
  t() %>% as.data.frame() %>%
  PROP() %>%
  INDEX() %>%
  t() %>% as.data.frame() %>%
  rownames_to_column(var = "dat_site") %>%
  pivot_longer(!dat_site, names_to = "LCT_shared", values_to = "index") %>%
  merge(e1,., by = c("dat_site", "LCT_shared"))%>%
  rename(read_count = count, read_index = index)%>%
  select(-c("date", "site", "hakai_link_date"))

#LCTs unique to each method ####
LCTs_eDNA <- unique(eDNA_long$LCT_shared)[!unique(eDNA_long$LCT_shared) %in% unique(bs_long$LCT_shared)]
LCTs_bs <- unique(bs_long$LCT_shared)[!unique(bs_long$LCT_shared) %in% unique(eDNA_long$LCT_shared)]

long <- merge(bs_long, eDNA_long, by = c("dat_site", "LCT_shared"), all.x = T, all.y = T) %>% replace(is.na(.), 0)

#were they detected in the same surveys
long$error <- if_else(long$abundance == 0 & long$read_index == 0, "not detected with either method",
                       if_else(long$abundance == 0, "eDNA detection only", 
                               if_else(long$read_index == 0, "beach seine detection only", "detected with both methods")))
#were they detected with either method at all
long <- long %>%
  mutate(gamma_error = if_else(LCT_shared %in% LCTs_eDNA, "eDNA only",
                               if_else(LCT_shared %in% LCTs_bs, "BS only", "both")))

a1 <- long %>%
  group_by(LCT_shared) %>%
  summarise(max_grams_per_m2 = max(grams_per_m2), sum_grams_per_m2 = sum(grams_per_m2))
a2 <- merge(long, a1, by = c("LCT_shared"))
long <- a2 %>%
  mutate(grams_per_m2_rel = grams_per_m2 / max_grams_per_m2)%>%
  mutate(grams_per_m2_perc = grams_per_m2 / sum_grams_per_m2)


write_rds(long, "Data/2022_10_31/derived_data/shared_long.rds")


#summaries of prevalence and abundance ####

# table for annotating ####
b1 <- distinct(long[c("LCT_shared", "gamma_error")]) %>%
  merge(., taxonomy_shared[c("LCT_shared", "all_species")], by = "LCT_shared")

b2 <- long %>%
  group_by(dat_site) %>%
  summarise(tot_biomass = sum(grams_per_m2)) %>%
  merge(long, ., by = "dat_site") %>%
  mutate(prop_biodens = grams_per_m2 / tot_biomass) %>%
  group_by(LCT_shared) %>%
  summarise(mean_biodensBS = mean(grams_per_m2), 
            max_biodensBS = max(grams_per_m2),
            total_biodensBS = sum(grams_per_m2),
            max_prop_biodensBS = max(prop_biodens),
            mean_read_eDNA = mean(read_count), 
            max_read_eDNA = max(read_count),
            prev_bs = sum(p_a_bs),
            prev_eDNA = sum(p_a_eDNA))

c1 <- long %>%
  filter(gamma_error == "both") %>%
  filter(error != "not detected with either method") %>%
  group_by(LCT_shared, error) %>%
  summarise(error_freq = length(error)) %>%
  pivot_wider(names_from = error, values_from = error_freq) %>%
  replace(is.na(.), 0) %>%
  `colnames<-` (c("LCT_shared","both","eDNA", "BS")) %>%
  mutate(both_per = both/(both + eDNA + BS),
         eDNA_per = eDNA/(both + eDNA + BS),
         BS_per = BS/(both + eDNA + BS))

c2 <- c1 %>%
  mutate(grp = case_when(both_per > 0.5 ~ "both_m",
                         eDNA_per > 0.5 ~ "eDNA_m",
                         BS_per > 0.5 ~ "BS_m")) %>%
  mutate(grp = ifelse(is.na(grp), "NP", grp)) %>%
  arrange(grp)

b3 <- merge(b1, b2, by = "LCT_shared")
b4 <- merge(b3, c2[c("LCT_shared", "grp")], by = "LCT_shared", all = T)

write.csv(b4,"./Figures_tables/LCT_table.csv")

# % total biomass of species not detected with BS
b5 <- b4 %>%
  group_by(gamma_error) %>%
  summarise(total_biodensBS = sum(total_biodensBS)) %>%
  mutate(percent_total_biomass = 100* total_biodensBS / sum(total_biodensBS))

#compare read index against abundance ####
long_shared <- long %>%
  filter(gamma_error == "both") %>%
  mutate(across(LCT_shared, factor, levels= c2$LCT_shared)) %>%
  merge(., c2[c("LCT_shared", "grp")], by = "LCT_shared")
  #  filter()
#All species
#Caption: Raw biomass against read index of species from paired beach seining and eDNA surveys (n = 46)
#from 12 location across 3 years on the Central Coast of British Columbia. Abundance is a sum across two
#beach seine set replicates for each site. Read index is a survey read abundance relative to the total 
#read abundance across that species allowing for a more direct comparison. Black points indicate surveys
#where the species was detected in both methods. Orange and red points indicate surveys where a species was
#only detected in beach seining or eDNA surveys, respectively.
#biomass
biomass_x_read_facet_spec <- 
  ggplot(long_shared, aes(y = weight, x = read_index, colour = error)) +
  geom_point(size = 0.9) +
  scale_color_manual(name = NULL,
                     values=c("detected with both methods" = "#000000", 
                              "beach seine detection only" = "#E69F00",
                              "eDNA detection only" = "#D55E00",
                              "not detected with either method" = "light grey")) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10), 
                     breaks = c(0,1,10,100,1000,10000,100000),
                     labels = c("0 ",
                                "1   ", 
                                expression(10^{"  "}), 
                                expression(10^{"2"}), 
                                expression(10^{"3"}), 
                                expression(10^{"4"}), 
                                expression(10^{"5"}))) +
  scale_x_continuous(breaks = c(0,1)) +
  labs(x = "Read Index (eDNA)",
       y = "Biomass (grams)") +
  #  theme(plot.title = element_text(size = 2)) +
  #  annotate("text", label = table(long$error)) +
  facet_wrap(~LCT_shared, ncol = 7, labeller = label_wrap_gen(width=1)) +
  theme_bw() +
  theme(legend.position = c(0.65,0.055),
        #      legend.spacing.y = unit(0, "cm"),
        #     legend.key.size = unit(0, 'cm'),
        #    legend.key.height = unit(0, "cm"),
        legend.direction="horizontal",
        legend.text=element_text(size=11),
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 11),
        strip.background = element_rect(fill = "white", colour = "gray", size = 1),
        #        margin(r = 1, unit = "pt"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text.x = element_text(size = 7)) +
  #  annotate("text", x = 0.5, y = 1, label = "Text No. 1")
  guides(colour = guide_legend(nrow=2, byrow=TRUE, override.aes = list(size=3)))

  biomass_x_read_facet_spec 

ggsave("./Figures_tables/biomass_x_read_facet_spec.png", 
       plot = biomass_x_read_facet_spec,
       width = 7.75, height = 7, units = "in")


#for biomass density 
long_shared <- long %>%
  filter(gamma_error == "both") %>%
  mutate(across(LCT_shared, factor, levels= c2$LCT_shared)) %>%
  merge(., c2[c("LCT_shared", "grp")], by = "LCT_shared")#%>%
#  filter(grams_per_m2 <10)

 
#linear models for individual species 
no_0s_shared <- long_shared %>%
  filter(error != "not detected with either method")

sig_LITs <- no_0s_shared %>%
  mutate(LCT_shared = as.character(LCT_shared)) %>%
  group_by(LCT_shared) %>%
  summarise(obs = length(LCT_shared)) %>%
  filter(obs > 15) 

LITs <- sig_LITs$LCT_shared
models <- list()

for (i in LITs) {
m1 <- filter(no_0s_shared, LCT_shared == i)
m2 <- lm(grams_per_m2_rel ~ read_index, data = m1) #change between relative and absulute biomass density here
m3 <- c(length(m1$LCT_shared), summary(m2)$coefficients[1, 1],summary(m2)$coefficients[2, 1], summary(m2)$r.squared)
models[[length(models)+1]] = m3
}

m4 <- as.data.frame(models)
colnames(m4) = LITs
m5 <- as.data.frame(t(m4))
colnames(m5) = c("n_observations", "intercept", "slope", "R2")
m6 <- rownames_to_column(m5, "LCT")
ggplot(m6, aes(x = R2)) +
  geom_histogram(binwidth = 0.05,colour = "black", fill = "grey90") +
  theme_classic() + 
  xlab(bquote(R^2)) +
  ylab("count (n=30)")

write_csv(m6, "./Figures_tables/BM_read_regressions.csv")

#mean and median of R2 vales for biomass density against abundance
mean(m6$R2) 
median(m6$R2)

#non-shared fishes ####

# identify fish undetected in each method
  long_bs_only <- long %>%
    filter(gamma_error == "BS only")
  long_eDNA_only <- long %>%
    filter(gamma_error == "eDNA only")
#creat order for plotting  
bs_unique_order <- long_bs_only %>%
  group_by(LCT_shared) %>%
  summarise(mean_abund = mean(abundance)) %>%
  .[order(.$mean_abund, decreasing = T),]
eDNA_unique_order <- long_eDNA_only %>%
  group_by(LCT_shared) %>%
  summarise(mean_abund = mean(read_index)) %>%
  .[order(.$mean_abund, decreasing = T),]

bs_only_abund <-
ggplot(long_bs_only, aes(x = LCT_shared, y = grams_per_m2)) +
  geom_jitter(shape=16, position=position_jitter(0.1), colour = "#E69F00") +
  scale_x_discrete(limits = bs_unique_order$LCT_shared) +
  theme_classic()+
  scale_y_continuous(trans = pseudo_log_trans(base = 10), 
                     breaks = c(0,1,10,100,1000,10000),
                     labels = c("0 ","1   ", expression(10^{"  "}), expression(10^{"2"}), expression(10^{"3"}), expression(10^{"4"}))) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2, size = 8.5),
        axis.title.y = element_text(size = 8)) +
  labs(x = NULL,
       y = "Grams / m^2")
bs_only_abund

eDNA_only_abund <-
ggplot(long_eDNA_only, aes(x = LCT_shared, y = read_index)) +
  geom_jitter(shape=16, position=position_jitter(0.1), colour = "#D55E00") +
  scale_x_discrete(limits = eDNA_unique_order$LCT_shared) +
  theme_classic() +
  labs(x = NULL,
       y = "Read Index") +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2, size = 8.5),
        axis.title.y = element_text(size = 8)) 
 # coord_flip()
eDNA_only_abund

#weight disribution ####
#plot weight distributions of from beach seining with colour coded observations
w7 <- long
w7$group1 <- factor(w7$error, levels = c("detected with both methods",
                                         "not detected with either method", #assign factors
                                          "eDNA detection only",
                                          "beach seine detection only"))
w7 <- with(w7, w7[order(group1),])                                          #assign order
#get counts of each observation type - add these to labels in plot below
sh1 <- w7 %>%
  group_by(group1) %>%
  summarise(count = length(group1))
#order taxa by mean biomass               
W_order <- w7 %>%
  group_by(LCT_shared) %>%
  summarise(mean_abund = mean(weight)) %>%
  .[order(.$mean_abund, decreasing = T),]
w8 <- w7 %>%
  mutate(weight1 = if_else(weight >= 0.8, weight, 
                              if_else(weight > 0, 0.8, 0)))
biomass_x_read_all <-
ggplot(w8, aes(x = LCT_shared, y = weight1, colour = group1)) +
  geom_jitter(shape=16, position=position_jitter(width = 0.3, height = 0.05), size = 2) +
  scale_color_manual(name = "Observation type", 
                     labels = c("detected with both methods (n = 386)",
                                "not detected with either method (n = 2517)",
                                "beach seine detection only (n = 233)",
                                "eDNA detection only (n = 682)"),
                     values=c("detected with both methods" = "#000000", 
                              "beach seine detection only" = "#E69F00",
                              "eDNA detection only" = "#D55E00",
                              "not detected with either method" = "lightgrey")) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10), 
                     breaks = c(0,1,10,100,1000,10000,100000),
                     labels = c("0","1   ", "10 ", expression(10^{"2"}), expression(10^{"3"}), expression(10^{"4"}), expression(10^{"5"})),
                     limits = c(-0.6, 100000)) +
  scale_x_discrete(limits = W_order$LCT_shared) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2, size = 8.5),
        axis.title.y = element_text(size = 14),
        legend.position = c(0.6,0.8),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(x = NULL,
       y = "Biomass (grams)") + 
#  theme(axis.text.x = element_text(face = W_order$group2)) +
  annotate("point", x = bs_unique_order$LCT_shared, y = -0.6,colour = "#E69F00", size = 2.5, alpha=1, shape = "triangle")+
  annotate("point", x = eDNA_unique_order$LCT_shared, y = -0.6, colour = "#D55E00", size = 2.5, alpha=1, shape = "triangle") +
  guides(colour = guide_legend(override.aes = list(size=3)))

  biomass_x_read_all
  
ggsave("./Figures_tables/biomass_x_read_all.png", 
       plot = biomass_x_read_all,
       width = 13.2, height = 7.5, units = "in")


#density disribution ####
#plot weight distributions of from beach seining with colour coded observations
q7 <- long
q7$group1 <- factor(q7$error, levels = c("detected with both methods", #assign factors,
                                         "not detected with either method",
                                         "eDNA detection only",
                                         "beach seine detection only"))
q7 <- with(q7, q7[order(group1),])                                          #assign order
#get counts of each observation type - add these to labels in plot below
sh1 <- q7 %>%
  group_by(group1) %>%
  summarise(count = length(group1))
#order taxa by mean biomass               
W_order <- q7 %>%
  group_by(LCT_shared) %>%
  summarise(mean_abund = mean(grams_per_m2), mean_eDNA_pa = sum(p_a_eDNA)) %>%
  .[order(.$mean_eDNA_pa, decreasing = T),] %>%
  .[order(.$mean_abund, decreasing = T),]

q9 <- q7 %>%
  filter(group1 == "detected with both methods" | group1 == "beach seine detection only" & grams_per_m2 > 0.00001)
q10 <- q7 %>%
  filter(group1 == "eDNA detection only" | group1 == "not detected with either method")
q11 <- q7 %>%
  filter(grams_per_m2 <0.00001 & grams_per_m2 > 0)

biomassDens_x_read_all <-
  ggplot() +
  geom_jitter(data = q9, aes(x = LCT_shared, y = grams_per_m2, colour = group1),shape=16, position=position_jitter(width = 0.25, height = 0.02), size = 1.5) +
  geom_jitter(data = q10, aes(x = LCT_shared, y = 0.00005, colour = group1),shape=16, position=position_jitter(width = 0.25, height = 0.1), size = 1.5) +
  geom_point(data = q11, aes(x = LCT_shared, y = 0.0001, colour = group1),shape=16, size = 1.5) +
  scale_color_manual(name = "Observation type", 
                     labels = c( "beach seine detection only (n = 233)",
                                 "detected with both methods (n = 386)",
                                "eDNA detection only (n = 682)",
                                "not detected with either method (n = 2517)"),
                     values=c("beach seine detection only" = "#E69F00",
                              "detected with both methods" = "#000000", 
                              "eDNA detection only" = "#D55E00",
                              "not detected with either method" = "lightgrey")) +
  scale_y_continuous(trans = 'log10',limits = c(0.000025,1e3), 
                     breaks = c(0.0001, 0.0001, 0.001,0.01,0.1,1,10,100,1000),
                     labels = c(expression(10^{"-5"}),
                                expression(10^{"-4"}),
                                expression(10^{"-3"}),
                                expression(10^{"-2"}),
                                expression(10^{"-1"}),
                                "1", 
                                "10",
                                expression(10^{"2"}), 
                                expression(10^{"3"}))) +
  scale_x_discrete(limits = W_order$LCT_shared) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2, size = 10),
        axis.title.y = element_text(size = 14),
        legend.position = c(0.3,0.8),
        legend.text=element_text(size = 12),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) +
  labs(x = NULL) +
 #      y =  "Biomass density (grams / m^2)") + 
  ylab(bquote('Biomass density '(g / m^2))) +
  #  theme(axis.text.x = element_text(face = W_order$group2)) +
  annotate("point", x = bs_unique_order$LCT_shared, y = 0.000025,colour = "#E69F00", size = 2.5, alpha=1, shape = "triangle")+
  annotate("point", x = eDNA_unique_order$LCT_shared, y = 0.000025, colour = "#D55E00", size = 2.5, alpha=1, shape = "triangle") +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  annotate("text", x = 14.2, y = 0.00023, label = expression(10^{"-9"}), size = 4) +
  annotate("segment", x= 13, y = 0.0003, xend = 13, yend = 0.00012, arrow = arrow(length = unit(.15,"cm"))) +
  annotate("text", x = 57.2, y = 0.00023, label = expression(10^{"-6"}), size = 4) +
  annotate("segment", x= 56, y = 0.0003, xend = 56, yend = 0.00012, arrow = arrow(length = unit(.15,"cm")))
biomassDens_x_read_all

ggsave("./Figures_tables/biomassDensity_x_read_all.png", 
       plot = biomassDens_x_read_all,
       width = 13.2, height = 7.5, units = "in")

biomassDens_x_read_all_small <-
  ggplot() +
  geom_jitter(data = q9, aes(x = LCT_shared, y = grams_per_m2, colour = group1),shape=16, position=position_jitter(width = 0.25, height = 0.02), size = 1.5) +
  geom_jitter(data = q10, aes(x = LCT_shared, y = 0.00005, colour = group1),shape=16, position=position_jitter(width = 0.25, height = 0.1), size = 1.5) +
  geom_point(data = q11, aes(x = LCT_shared, y = 0.0001, colour = group1),shape=16, size = 1.5) +
 # geom_point(data = bs_unique_order, aes(x = LCT_shared, y = 0.000025),colour = "#E69F00", size = 2, alpha=1, shape = "triangle")+
  scale_color_manual(name = "Observation type", 
                     labels = c( "beach seine detection only (n = 233)",
                                 "detected with both methods (n = 386)",
                                 "eDNA detection only (n = 682)",
                                 "not detected with either method (n = 2517)"),
                     values=c("beach seine detection only" = "#E69F00",
                              "detected with both methods" = "#000000", 
                              "eDNA detection only" = "#D55E00",
                              "not detected with either method" = "lightgrey")) +
  scale_y_continuous(trans = 'log10',limits = c(0.000005,1e3),
                     expand = c(0, 0),
                     breaks = c(0.0001, 0.001,0.01,0.1,1,10,100,1000),
                     labels = c(expression("<" ~10^{"-4"}),
                                expression(10^{"-3"}),
                                expression(10^{"-2"}),
                                expression(10^{"-1"}),
                                "1", 
                                "10",
                                expression(10^{"2"}), 
                                expression(10^{"3"}))) +
  scale_x_discrete(limits = W_order$LCT_shared) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 11),
        legend.position = c(0.75,0.8),
        legend.text=element_text(size = 12),
        plot.margin = margin(0.5, 0.5, 2, 0.5, "cm")) +
  labs(x = "taxa in ranked order of mean biomass") +
  #      y =  "Biomass density (grams / m^2)") + 
  ylab(bquote('biomass density '(g / m^2))) +
  #  theme(axis.text.x = element_text(face = W_order$group2)) +
  annotate("point", x = bs_unique_order$LCT_shared, y = 0.000025,colour = "#E69F00", size = 2, alpha=1, shape = "triangle")+
  annotate("point", x = eDNA_unique_order$LCT_shared, y = 0.000025, colour = "#D55E00", size = 2, alpha=1, shape = "triangle") +
  guides(colour = guide_legend(override.aes = list(size=3))) +
 # annotate("text", x = 14.8, y = 0.00023, label = expression(10^{"-9"}), size = 4) +
 # annotate("segment", x= 13, y = 0.0003, xend = 13, yend = 0.00012, arrow = arrow(length = unit(.15,"cm"))) +
 # annotate("text", x = 57.8, y = 0.00023, label = expression(10^{"-6"}), size = 4) +
  #annotate("segment", x= 56, y = 0.0003, xend = 56, yend = 0.00012, arrow = arrow(length = unit(.15,"cm"))) +
  geom_hline(yintercept = 0.0001, colour = "gray", linetype = "dotted")#+
  #annotate("text", x = 20, y = 0.000005, label = "text")  +
#  coord_cartesian(ylim = c(0.00002, 1e3), clip = "off")

biomassDens_x_read_all_small

ggsave("./Figures_tables/biomassDensity_x_read_all_small.png", 
       plot = biomassDens_x_read_all_small,
       width = 8, height = 4, units = "in")

#euler plot
library(eulerr)
#eDNA 
628 +386
#1014
#BS
233+386
#619
#both
#396
fit <- euler(c("A" = 1014, "B" = 619, "A&B" = 396))
plot(fit)




#summaries of prevalence and abundance

# table for annotating ####
b1 <- distinct(long[c("LCT_shared", "gamma_error")]) %>%
  merge(., taxonomy_shared[c("LCT_shared", "all_species")], by = "LCT_shared")

b2 <- long %>%
  group_by(dat_site) %>%
  summarise(tot_biomass = sum(grams_per_m2)) %>%
  merge(long, ., by = "dat_site") %>%
  mutate(prop_biodens = grams_per_m2 / tot_biomass) %>%
  group_by(LCT_shared) %>%
  summarise(mean_biodensBS = mean(grams_per_m2), 
            max_biodensBS = max(grams_per_m2), 
            max_prop_biodensBS = max(prop_biodens),
            mean_read_eDNA = mean(read_count), 
            max_read_eDNA = max(read_count),
            prev_bs = sum(p_a_bs),
            prev_eDNA = sum(p_a_eDNA))

c1 <- long %>%
  filter(gamma_error == "both") %>%
  filter(error != "not detected with either method") %>%
  group_by(LCT_shared, error) %>%
  summarise(error_freq = length(error)) %>%
  pivot_wider(names_from = error, values_from = error_freq) %>%
  replace(is.na(.), 0) %>%
  `colnames<-` (c("LCT_shared","both","eDNA", "BS")) %>%
  mutate(both_per = both/(both + eDNA + BS),
         eDNA_per = eDNA/(both + eDNA + BS),
         BS_per = BS/(both + eDNA + BS))

c2 <- c1 %>%
  mutate(grp = case_when(both_per > 0.5 ~ "both_m",
                       eDNA_per > 0.5 ~ "eDNA_m",
                       BS_per > 0.5 ~ "BS_m")) %>%
  mutate(grp = ifelse(is.na(grp), "NP", grp))

b3 <- merge(b1, b2, by = "LCT_shared")
b4 <- merge(b3, c2[c("LCT_shared", "grp")], by = "LCT_shared", all = T)

write.csv(b4,"./Figures_tables/LCT_table.csv")



#summary of abundant fishes: how many were missed by eDNA when biomass density was > 1 g/m^2
o1 <- q7 %>%
  filter(grams_per_m2 > 1) %>%
  group_by(error) %>%
  summarise(dd = length((error)))

# calculate % above 1 g/m^2
s6 <- q7 %>%
  filter(error == "detected with both methods" | error == "beach seine detection only") #%>%
  #filter(grams_per_m2 > 1)

  
# bin by biomass for plot following

l1 <- q7 %>%
  filter(p_a_bs != 0) %>%
  dplyr::select(c("p_a_bs", "p_a_eDNA", "grams_per_m2")) %>%
  mutate(bin = cut_number(grams_per_m2, n=10)) %>%
  mutate(bin = as.numeric(factor(bin))) 

l2 <- l1 %>%
  group_by(bin, p_a_eDNA) %>%
  summarise(g = length(p_a_eDNA)) %>%
  pivot_wider(names_from = p_a_eDNA, values_from = g) %>%
  `colnames<-`(c("bin", "eDNA_a", "eDNA_p")) %>%
  mutate(agree_perc = eDNA_p / (eDNA_a + eDNA_p) *100)

l3 <- l1 %>%
  group_by(bin) %>%
  summarise(mean = mean(grams_per_m2),
            min = min(grams_per_m2),
            max = max(grams_per_m2))

l4 <- merge(l2, l3, by = "bin")

missed_obs_ebe <-
  ggplot(l4, aes(x = mean, y = agree_perc)) +
  geom_point() +
  #  geom_text(vjust = 0.1, nudge_y = 1.5) +
  theme_classic()+ 
  ylim(0,100) +
  labs(y = "% of positive beach seine \n detections detected with eDNA") +
  xlab(bquote('Biomass density '(g / m^2)))   +
  scale_x_continuous(trans = pseudo_log_trans(base = 10), 
                     breaks = c(0,0.1,1,10,100),
                     limits=c(0, 270),
                     labels = c("0       ", "0.1", "1.0", "10", '100')) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  geom_vline(xintercept = c(0, 0.1, 1, 10), colour = "gray",linetype = "dashed") +
  geom_pointrange(aes(xmax = max, xmin = min))

missed_obs_ebe
ggsave("./Figures_tables/missed_obs_ebe.png", 
       plot = missed_obs_ebe,
       width = 6, height = 3, units = "in")


#SADs

#envonmental data
env <- read_rds("Data/2022_10_31/derived_data/environmental.rds")

long1 <- long %>%
  merge(., env[c("dat_site","h100m", "silt_percent", "site", "year", "month")], by = "dat_site")

t1 <- long1 %>%
  group_by(dat_site) %>%
  mutate(rank_bio = order(order(grams_per_m2, decreasing=TRUE)),
         rank_abund = order(order(abundance, decreasing=TRUE))) %>%
  filter(grams_per_m2 > 0.00000000000000000000001) %>%
  mutate(p_a_eDNA = as.character(p_a_eDNA)) %>%
  mutate(grams_per_m2 = ifelse(grams_per_m2 <= 0.0001, 0.0001, grams_per_m2))

t2 <- t1 %>% filter(p_a_eDNA == 0)

t3 <- t1 %>% filter(p_a_eDNA == 1)
t4 <- l4 %>% mutate(agree_perc = round(agree_perc),
                    min = ifelse(min <=0.0001, 0.0001, min))



pd <- position_dodge(0)
ggplot() +
  #  geom_point(size = 2) +
  geom_line(data = t1, aes(x = rank_bio, y = grams_per_m2, group = dat_site, colour = h100m),position = pd)+
  # scale_colour_gradient(name = "habitat richness\n100-m radius", low = "darkgreen", high = "yellow")+
  scale_colour_gradientn(name = "habitat richness\n100-m radius", colours =  c("green", "orange", "red"))+
  
  geom_point(data = t2, aes(x = rank_bio, y = grams_per_m2, shape = p_a_eDNA),position = pd)+
  scale_shape_manual(name = "eDNA false\nnegative", labels =c(""), values=c(1, 3))+
  scale_y_continuous(trans = 'log10',limits = c(0.0001,1e3), 
                     breaks = c(0.001,0.01,0.1,1,10,100,1000),
                     labels = c(expression(10^{"-3"}),
                                expression(10^{"-2"}),
                                expression(10^{"-1"}),
                                "1", 
                                "10",
                                expression(10^{"2"}), 
                                expression(10^{"3"}))) +
  theme_classic()+ 
  theme(legend.position = c(0.8,0.7))  + 
  ylab(bquote('Biomass density '(g / m^2))) +
  xlab("Rank")

pd <- position_dodge(0)
ggplot() +
  #  geom_point(size = 2) +
  geom_line(data = t1, aes(x = rank_bio, y = grams_per_m2, group = dat_site, colour = p_a_eDNA),position = pd, colour = "lightgrey")+
 # scale_colour_gradient(name = "habitat richness\n100-m radius", low = "darkgreen", high = "yellow")+
#  scale_colour_gradientn(name = "habitat richness\n100-m radius", colours =  c("green", "orange", "red"))+
  geom_jitter(data = t3, aes(x = rank_bio, y = grams_per_m2, shape = p_a_eDNA), size = 1.5, width = 0.05, height = 0, colour = "black", alpha = 0.3)+
  geom_jitter(data = t2, aes(x = rank_bio, y = grams_per_m2, shape = p_a_eDNA), size = 2, width = 0.05, height = 0, colour = "#D55E00", alpha = 0.3)+
  
   scale_shape_manual(name = "eDNA\nobservation\ntype", labels =c("false negative", "true positive"), values=c(19, 20))+
  geom_hline(yintercept = 0.0001, colour = "gray", linetype = "dotted") +
  scale_y_continuous(trans = 'log10',limits = c(0.0001,1e3), 
                     breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000),
                     labels = c(expression("<" ~10^{"-4"}),
                                expression(10^{"-3"}),
                                expression(10^{"-2"}),
                                expression(10^{"-1"}),
                                "1", 
                                "10",
                                expression(10^{"2"}), 
                                expression(10^{"3"}))) +
  theme_classic()+ 
  theme(legend.position = c(0.7,0.6))  + 
  ylab(bquote('Biomass density '(g / m^2))) +
  xlab("Rank") +
  geom_ysidedensity(data = t1, aes(y = grams_per_m2, fill = p_a_eDNA), alpha = 0.8) + 
  scale_ysidex_continuous(breaks = NULL, labels = NULL) +
  geom_xsidedensity(data = t1, aes(x = rank_bio, fill = p_a_eDNA), alpha = 0.8) + 
  scale_xsidey_continuous(breaks = NULL, labels = NULL) +
  scale_fill_manual(name = "dd",values=c("#D55E00", "lightgrey"))
  



pd <- position_dodge(0)
ggplot() +
  #  geom_point(size = 2) +
  geom_line(data = t1, aes(x = rank_abund, y = count_per_m2, group = dat_site, colour = p_a_eDNA),position = pd, colour = "lightgrey")+
  # scale_colour_gradient(name = "habitat richness\n100-m radius", low = "darkgreen", high = "yellow")+
  #  scale_colour_gradientn(name = "habitat richness\n100-m radius", colours =  c("green", "orange", "red"))+
  geom_jitter(data = t3, aes(x = rank_abund, y = count_per_m2, shape = p_a_eDNA), size = 1.5, width = 0.05, height = 0, colour = "black", alpha = 0.3)+
  geom_jitter(data = t2, aes(x = rank_abund, y = count_per_m2, shape = p_a_eDNA), size = 2, width = 0.05, height = 0, colour = "#D55E00", alpha = 0.3)+
  
  scale_shape_manual(name = "eDNA\nobservation\ntype", labels =c("false negative", "true positive"), values=c(19, 20))+
  geom_hline(yintercept = 0.0001, colour = "gray", linetype = "dotted") +
  scale_y_continuous(trans = 'log10',limits = c(0.0001,1e3), 
                     breaks = c(0.0001,0.001,0.01,0.1,1,10,100,1000),
                     labels = c(expression("<" ~10^{"-4"}),
                                expression(10^{"-3"}),
                                expression(10^{"-2"}),
                                expression(10^{"-1"}),
                                "1", 
                                "10",
                                expression(10^{"2"}), 
                                expression(10^{"3"}))) +
  theme_classic()+ 
  theme(legend.position = c(0.7,0.6))  + 
  ylab(bquote('Biomass density '(g / m^2))) +
  xlab("Rank") +
  geom_ysidedensity(data = t1, aes(y = count_per_m2, fill = p_a_eDNA), alpha = 0.8) + 
  scale_ysidex_continuous(breaks = NULL, labels = NULL) +
  geom_xsidedensity(data = t1, aes(x = rank_abund, fill = p_a_eDNA), alpha = 0.8) + 
  scale_xsidey_continuous(breaks = NULL, labels = NULL) +
  scale_fill_manual(name = "dd",values=c("#D55E00", "lightgrey"))

  
  
  

  missed_obs_ebe2 <-
    ggplot(t4, aes(y = mean, x = agree_perc)) +
    geom_point() +
    #  geom_text(vjust = 0.1, nudge_y = 1.5) +
    theme_classic()+ 
    xlim(25,100) +
    labs(x = "% of beach seine detections\ndetected with eDNA") +
    ylab(bquote('Biomass density '(g / m^2)))   +
    scale_y_continuous(trans = 'log10',limits = c(0.0001,1e3), 
                       breaks = c(0.001,0.01,0.1,1,10,100,1000),
                       labels = c(expression(10^{"-3"}),
                                  expression(10^{"-2"}),
                                  expression(10^{"-1"}),
                                  "1", 
                                  "10",
                                  expression(10^{"2"}), 
                                  expression(10^{"3"}))) +
    theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=0.5)) +
  #  geom_hline(yintercept = c(0, 0.1, 1, 10), colour = "gray",linetype = "dashed") +
    geom_pointrange(aes(ymax = max, ymin = min))
  
  missed_obs_ebe2


x1_1 <- t1 %>% group_by(dat_site) %>% summarise(richness = length(dat_site)) 
  
x1 <- t1 %>%
  group_by(dat_site, error, h100m, ) %>%
  summarise(count = length(dat_site)) %>%
  pivot_wider(names_from = error, values_from = count)  %>%
  `colnames<-`(c("dat_site", "h100m","bs","both")) %>%
  mutate(bs = if_else(is.na(bs), 0, bs)) %>%
  mutate(perc_missed = bs / (bs + both) *100,
         bs_rich = bs + both) 

write_rds(x1, "Data/2022_10_31/percent_missed.rds")

x2 <- x1 %>%
  merge(., env[c("dat_site", "site")], by = "dat_site") %>%
  group_by(site, h100m) %>%
  summarise(bs_mean = mean(bs),
            both_mean = mean(both))%>%
  mutate(perc_missed = bs_mean / (bs_mean + both_mean) *100,
         bs_rich = bs_mean + both_mean) 

x3 <- x1 %>%
  merge(., env[c("dat_site", "site")], by = "dat_site") %>%
  group_by(site, h100m) %>%
  summarise(perc_missed = mean(perc_missed))

ggplot(x1, aes(x = bs_rich, y = perc_missed)) +
  geom_point() +
  theme_classic() +
  xlim(1,23) +
  xlab("beach seine richness") +
  ylab("% missed\nby eDNA") 
  


ggplot(x2, aes(x = bs_rich, y = perc_missed)) +
  geom_jitter(height = 0.01) +
  theme_classic() +
  xlim(0,23) +
  xlab("Beach seine richness") +
  ylab("% missed\nby eDNA")

ggplot(x1, aes(x = h100m, y = perc_missed)) +
  geom_jitter(height = 0.01) +
  theme_classic() +
  xlim(0,4) +
  xlab("Habitat richness") +
  ylab("% missed\nby eDNA")


ggplot(x3, aes(x = h100m, y = perc_missed)) +
  geom_jitter(height = 0.01) +
  theme_classic() +
  xlim(0,4) +
  xlab("beach seine richness") +
  ylab("% missed\nby eDNA")


#percent missed by eDNA correlation with beach seine richness
cor.test(x = x2$bs_rich, y= x2$perc_missed, method = "spearman")

#plot richness in each method against richness at both 1000m and 100m ####

#calculate richness 
rich <- long %>%
  dplyr::select(c("dat_site", "p_a_bs", "p_a_eDNA"))  %>%
  group_by(dat_site) %>%
  mutate(across(c("p_a_bs", "p_a_eDNA"), sum)) %>%
  `colnames<-` (c("dat_site","rich_b","rich_e")) %>%
  distinct() %>%
  mutate(rich_diff = rich_e - rich_b) %>%
  merge(., env[c("dat_site","h100m", "h1000m", "silt_percent")], by = "dat_site")

rich1 <- rich %>% 
  pivot_longer(cols = c("h100m", "h1000m", "silt_percent"), names_to = "scale", values_to = "hab_rich")

rich2 <- rich %>%
  pivot_longer(cols = c("rich_e", "rich_b"), names_to = "method", values_to = "richness")%>%
  pivot_longer(cols = c("h100m", "h1000m", "silt_percent"), names_to = "scale", values_to = "hab_rich") %>% 
  mutate(site = gsub( " .*$", "", dat_site))

h100m <- rich2 %>% filter(scale == "h100m")
h1000m <- rich2 %>% filter(scale == "h1000m")
silt_percent <- rich2 %>% filter(scale == "silt_percent") %>%
  mutate(hab_rich = hab_rich * 100)

h100m_e <- rich2 %>% filter(scale == "h100m" & method == "rich_e")
h100m_b <- rich2 %>% filter(scale == "h100m" & method == "rich_b")
h1000m_e <- rich2 %>% filter(scale == "h1000m" & method == "rich_e")
h1000m_b <- rich2 %>% filter(scale == "h1000m" & method == "rich_b")

silt_percent_e <- rich2 %>% filter(scale == "silt_percent" & method == "rich_e") %>%
  mutate(hab_rich = hab_rich * 100)
silt_percent_b <- rich2 %>% filter(scale == "silt_percent" & method == "rich_b") %>%
  mutate(hab_rich = hab_rich * 100)


rich_b_h100m <- glm(richness ~ hab_rich, data = h100m_b)
summary(rich_b_h100m)


p_h100m <- ggplot() +
  geom_jitter(data = h100m_e, aes(x = hab_rich, y = richness), shape = 16, width = 0.2, colour = "#D55E00", size = 2.5)+
  geom_smooth(data = h100m_e, aes(x = hab_rich, y = richness),method = lm, se = T, colour = "#D55E00")+
  geom_jitter(data = h100m_b, aes(x = hab_rich, y = richness), shape = 17, width = 0.2, colour = "#E69F00", size = 2.5)+
  geom_smooth(data = h100m_b, aes(x = hab_rich, y = richness),method = lm, se = T, colour = "#E69F00")+
#  geom_line(data = h100m, aes(x = hab_rich, y = richness, group = site)) +
  theme_classic() +
  xlab("habitat richness (100-m radius)") +
  ylab("taxonomic richness") +
  annotate("text", x = -0.1, y = 37, label = "B", fontface =2, size = 5)
p_h100m

p_h1000m <- ggplot() +
  geom_jitter(data = h1000m_e, aes(x = hab_rich, y = richness), shape = 16, width = 0.2, colour = "#D55E00", size = 2.5)+
  geom_smooth(data = h1000m_e, aes(x = hab_rich, y = richness),method = lm, se = T, colour = "#D55E00")+
  geom_jitter(data = h1000m_b, aes(x = hab_rich, y = richness), shape = 17, width = 0.2, colour = "#E69F00", size = 2.5)+
  geom_smooth(data = h1000m_b, aes(x = hab_rich, y = richness),method = lm, se = T, colour = "#E69F00")+
  #  geom_line(data = h100m, aes(x = hab_rich, y = richness, group = site)) +
  theme_classic() +
  xlim(2.5,5.5)+
  xlab("habitat richness (1000-m radius)") +
  ylab("taxonomic richness") +
  annotate("text", x = 0, y = 37, label = "B", fontface =2, size = 5)
p_h1000m

reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

p_silt <- ggplot() +
  geom_jitter(data = silt_percent_e, aes(x = hab_rich, y = richness), shape = 16, width = 0.01, colour = "#D55E00", size = 2.5)+
  geom_smooth(data = silt_percent_e, aes(x = hab_rich, y = richness),method = lm, se = T, colour = "#D55E00")+
  geom_jitter(data = silt_percent_b, aes(x = hab_rich, y = richness), shape = 17, width = 0.01, colour = "#E69F00", size = 2.5)+
  geom_smooth(data = silt_percent_b, aes(x = hab_rich, y = richness),method = lm, se = T, colour = "#E69F00")+
  theme_classic()+ 
  scale_x_continuous(trans=reverselog_trans(10)) +
 # scale_x_reverse() +
  xlab("% fine sediment (reversed log10 scale)") +
  ylab("taxonomic richness")+
  scale_shape_manual(values = c(17, 16)) +
  annotate("segment", x= 32,xend = 1,y=5, yend = 5, linewidth = 0.5, arrow = arrow(length = unit(0.3, 'cm'))) +
  annotate("text", x = 6, y = 6.2, label = "increased seawater movement") +
  annotate("pointrange", x= 4,xmin = 3.5, xmax = 4.5, y=34, linewidth = 1, colour = "#D55E00") +
  annotate("text", x = 8, y = 34, label = "eDNA")+
  annotate("pointrange", x= 4,xmin = 3.5, xmax = 4.5, y=32, linewidth = 1, colour = "#E69F00") +
  annotate("text", x = 11, y = 32, label = "beach seine") +
  annotate("text", x = 30, y = 37, label = "A", fontface =2, size = 5)
p_silt

plot_all <- 
  ggarrange(p_silt, p_h100m,  
            ncol = 2, nrow = 1,
            heights = c(1,1)) 
plot_all

ggsave("./Figures_tables/biomassDensity_x_read_all.png", 
       plot = plot_all,
       width = 8, height =4, units = "in")









