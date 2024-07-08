#Correct for mis-ids in beach seine surveys (clustering species into LITs (in code "LCT") for indistinguishable species.
#Resolve taxonomy for eDNA and beach seining - so both share same LITs

#packages
library(tidyverse)
library(readxl)
library(lubridate)
#read in raw data
bs <- read_csv("no_occupancy_model/Data/2022_10_31/bs_data/bs_abundance.csv")
taxa_BS <- readRDS("no_occupancy_model/Data/2022_10_31/derived_data/taxa_BS_resolved.rds")
env <- read_csv("no_occupancy_model/Data/2022_10_31/bs_data/bs_habitat.csv")
sitesurvey_data_eDNA <- readRDS("no_occupancy_model/Data/2022_10_31/derived_data/sitesurvey_data.rds")
taxa_eDNA <- read_rds("no_occupancy_model/Data/2022_10_31/taxonomy_eDNA.rds")

a1 <- bs %>%  #select surveys with paired eDNA
  mutate(date = ymd(paste(year, month, day, sep= ' '))) %>%
  filter(date %in% sitesurvey_data_eDNA$hakai_link_date)
a2 <- taxa_BS %>% #select taxa observed in paired surveys 
  filter(query %in% a1$species)
obs_to_remove <- filter(a2, is.na(genus) | query == "unro") # these taxa have very poor resolution and are very infrequent - remove from dataset
groups_BS <- filter(a2, is.na(species) & !is.na(genus) & query != "unro") %>% #taxa groups - these species were frequently misidentified in the field
  mutate(all_species = c("Sebastes melanops, Sebastes flavidus", 
                         "Sebastes caurinus, Sebastes maliger", 
                         "Artedius fenestralis, Artedius harringtoni, Artedius lateralis")) %>%
  mutate(LCT = c("Sebastes2", 
                 "Sebastes1", 
                 "Artedius1")) %>%  #assign groups used in eDNA or assign new groups %>%
  mutate(species = LCT) %>%
  mutate(genus = LCT) %>%
  mutate(level = genus) %>%
  separate(all_species, c("A", "B", "C", "D"), sep = ", ", remove = F) %>%
  pivot_longer(cols = c("A", "B", "C", "D")) %>%
  rename('spec' = 'value') %>%
  select(-c("name")) %>%
  filter(!is.na(spec)) %>%
  distinct() %>%
  .[c("spec", "level", "LCT", "class", "order", "family", "genus", "species", "all_species", "query")]
  
#groups from BS
BS_taxonomy <- a2[c(13, 3:6, 16)] %>%
  filter(!is.na(species)) %>%
  merge(groups_BS[c("spec", "LCT", "all_species")], ., by.x = "spec", by.y = "species", all.y = T) %>%
  mutate(LCT = if_else(is.na(LCT), spec, LCT)) %>%
  mutate(all_species = if_else(is.na(all_species), spec, all_species))
write_rds(BS_taxonomy, "Data/2022_10_31/taxonomy_BS.rds")

#select highest grouping for each taxa from either BS or eDNA
shared_LCT <- merge(BS_taxonomy[c("spec", "LCT")], taxa_eDNA[c("spec", "LCT")], by = "spec", all.x = T, all.y = T) %>%
  mutate(LCT_shared = if_else(is.na(LCT.x),LCT.y,
                              if_else(is.na(LCT.y),LCT.x,
                                      if_else(LCT.x == LCT.y, LCT.x, 
                                              if_else(grepl("\\d", LCT.x), LCT.x,
                                                      if_else(grepl("\\d", LCT.y), LCT.y, "XXXXXXXXX")))))) %>%
  mutate(LCT_shared = if_else(LCT_shared == "Ammodytes hexapterus" |       #manual edits for taxa we are confident are the same thing
                              LCT_shared == "Ammodytes personatus", 
                                "Ammodytes1", LCT_shared)) %>%
  mutate(LCT_shared = if_else(LCT_shared == "Citharichthys sordidus" | 
                                LCT_shared == "Citharichthys stigmaeus", 
                              "Citharichthys1", LCT_shared))


BS_taxonomy_reconciled <- merge(BS_taxonomy, shared_LCT[c("spec", "LCT_shared")], by = "spec") %>%
  rename(LCT_BS = LCT) %>%
  mutate(all_species = if_else(LCT_shared == "Citharichthys1", 
                              "Citharichthys sordidus, Citharichthys stigmaeus", all_species))%>%
  mutate(all_species = if_else(LCT_shared == "Ammodytes1", 
                               "Ammodytes hexapterus, Ammodytes personatus", all_species))

write_rds(BS_taxonomy_reconciled, "no_occupancy_model/Data/2022_10_31/derived_data/BS_taxonomy_reconciled.rds")

ASV_taxonomy <- read_rds("no_occupancy_model/Data/2022_10_31/derived_data/ASV_taxonomy_eDNA20221103.rds") 
ASV_taxonomy_reconciled <- ASV_taxonomy %>%
  separate(all_species, c("A", "B", "C", "D"), sep = ", ", remove = F) %>%
  pivot_longer(cols = c("A", "B", "C", "D")) %>%
  rename('spec' = 'value') %>%
  filter(!is.na(spec)) %>%
  merge(., shared_LCT[c("spec", "LCT_shared")], by = "spec") %>%
  rename(LCT_eDNA = LCT) %>%
  select(-c("name", "spec")) %>%
  mutate(species_shared = LCT_shared) %>%
  distinct()%>%
  mutate(all_species = if_else(LCT_shared == "Citharichthys1", 
                               "Citharichthys sordidus, Citharichthys stigmaeus", all_species))%>%
  mutate(all_species = if_else(LCT_shared == "Ammodytes1", 
                               "Ammodytes hexapterus, Ammodytes personatus", all_species))%>%
  mutate(all_species = if_else(LCT_shared == "Sebastes2", 
                               "Sebastes - alutianus, alutus, elongatus, flavidus, melanistictus, melanops,  ruberrimus", all_species))


write_rds(ASV_taxonomy_reconciled, "no_occupancy_model/Data/2022_10_31/derived_data/ASV_taxonomy_reconciled.rds")

taxa <- BS_taxonomy_reconciled

bs <- a1 %>% rename("query" = "species") #rename var

t1 <- merge(bs, taxa, by = "query") #merge species data with taxonomy

#make wide data (maintaining replicates) - to be used for occupancy modelling
t2 <- t1 %>%                                              #use this later for occupancy modelling (retains replicates)
  group_by(site, year, month, day, replicate, LCT_BS) %>%
  summarise(abundance = sum(abundance)) %>%
  spread(LCT_BS, abundance)

saveRDS(t2, "no_occupancy_model/Data/2022_10_31/derived_data/bs_specmatrix_retainreplicates.rds") # maybe we do some occupancy modelling later

#make wide data
t3 <- t1 %>%
  group_by(site, year, month, day, LCT_BS) %>%
  summarise(abundance = sum(abundance)) %>%
  spread(LCT_BS, abundance)

#create survey ID - I think this will improve filtering at later steps
t3$loc_dat_depth <- as.numeric(factor(paste0(t3$site, t3$year, t3$month, t3$day)))
t4 <- t3 %>%
  mutate(surveyID = paste("bss", loc_dat_depth, sep = "")) %>%
  relocate(surveyID, .after = day) 
t4$loc_dat_depth <- NULL

t5 <- t4[c(5:ncol(t4))] %>% #select species matrix with survey ID
  column_to_rownames("surveyID") 

t5[is.na(t5)] <- 0

saveRDS(t5, "no_occupancy_model/Data/2022_10_31/derived_data/bs_specmatrix.rds")

t6 <- t4[c("site", "year", "month", "day", "surveyID")] # select survey data
saveRDS(t6, "no_occupancy_model/Data/2022_10_31/derived_data/bs_surveydata.rds") #could add environmental to this at some point


#reconcile all_species for each method
s1 <- rbind(ASV_taxonomy_reconciled[c("LCT_shared", "all_species", "class", "order", "family", "genus")],
            BS_taxonomy_reconciled[c("LCT_shared", "all_species", "class", "order", "family", "genus")]) %>%
  distinct()
#identify LCTs that need amendments 
s2 <- s1 %>% 
  group_by(LCT_shared) %>%
  summarise(LCTs = length(LCT_shared)) %>%
  filter(LCTs > 1)

s3 <- filter(s1, LCT_shared %in% s2$LCT_shared) %>%
  mutate(length = nchar(all_species)) %>%
  group_by(LCT_shared) %>%
  filter(length == max(length)) %>%
  select(-c("length"))

s4 <- s1 %>%
  filter(!LCT_shared %in% s2$LCT_shared) %>%
  rbind(s3)

write_rds(s4, "no_occupancy_model/Data/2022_10_31/derived_data/taxonomy_combined.rds")


