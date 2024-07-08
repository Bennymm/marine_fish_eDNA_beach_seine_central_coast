

#packages
library(rfishbase)
library(tidyverse)
library(here)
library(taxize)
library(readxl)
library(lubridate)

#data relating total length to fork length
tr <- read_csv("no_occupancy_model/Data/2022_10_31/specandgen_traits.csv")
data <- filter(tr, species == "y" & replicate == 1) %>%
  mutate(id = sub("_", " ", id)) %>%
  rename(spec = species) %>%
  rename(species = id) %>%
  mutate(species = tolower(species))


return1 <- tol_resolve(data$species) 
return <- return1 %>%  
  mutate(species = if_else(unique_name == "Oncorhynchus mykiss (species in domain Eukaryota)", 
                           "Oncorhynchus mykiss", 
                           unique_name),
         species = if_else(unique_name == "Rhacochilus vacca", 
                           "Phanerodon vacca", 
                           unique_name)) %>%
  #remove species that have been merged into another
  filter(search_string != "bothrocara remigerum" &
         search_string != "inopsetta ischyra" &
         search_string != "lycodapus grossidens")

traits <- merge(data, return, by.x = "species", by.y = "search_string") %>%
  dplyr::select(!species) %>%
  rename(species = species.y)%>%
  relocate(species)

# calculate total length : fork length rations from trait measurements taken from book ####

m1 <- estimate(return$species)
m2 <- m1[c("Species", "a", "b")] %>%
  rename(species = Species)

#read in BS length and abundance data and taxize species names
bs_length <- read.csv("no_occupancy_model/Data/2022_10_31/bs_data/bs_lengths.csv") %>%
  rename("query" = "species") %>%
  mutate(query = tolower(query)) %>%
  mutate(site = tolower(site))
bs_abund <- read.csv("no_occupancy_model/Data/2022_10_31/bs_data/bs_abundance.csv") %>%
  rename(query = species) %>%
  group_by(site, year, month, day, query) %>%
  summarise(count = sum(abundance, na.rm = TRUE))%>%
  mutate(site = tolower(site))

#bs_codes <- read_rds("no_occupancy_model/Data/2022_10_31/derived_data/taxa_bs_resolved.rds")
  #replace query with closest relative (only for weight estimates) and write
#write_csv(bs_codes, "no_occupancy_model/Data/2022_10_31/derived_data/names_for_length.csv")

bs_codes <- read.csv("no_occupancy_model/Data/2022_10_31/derived_data/names_for_length_1.csv") %>%
  dplyr::select(c("species", "query")) %>%
  mutate(species = if_else(species == "Rhacochilus vacca", 
                    "Phanerodon vacca", 
                    species))


gt1 <- merge(bs_length, bs_codes, by = "query", all.x = TRUE, all.y = FALSE) #%>%
 # filter(species != "Cottus aleuticus")
gt2 <- merge(gt1, traits[c("species", "fl_tl")], by = "species", all.x = TRUE, all.y = FALSE)
gt2[c("comments", "entered", "checked")] <- NULL


gt3 <- merge(gt2, m2, by = "species", all.x = TRUE, all.y = FALSE)

gt4 <- gt3 %>%
  mutate(weight = a * (length / 10 * fl_tl) ^ b)  # W = a * L^b, L = length(mm) /10 * fl:tl ratio 
#nas <- filter(gt4, is.na(weight))

sp1 <- gt4[c("species", "query", "site", "year", "month", "day", "weight")] 
sp1$num <- as.character(as.numeric(factor(paste0(sp1$species, sp1$query, sp1$site, sp1$year, sp1$month, sp1$day))))
sp1$link = paste("lnk", sp1$num, sep = "")

#split into 2 groups: 1) all fish were measured, 2) a subset were measured, or maybe 3) none measured
sp2 <- sp1 %>%
  group_by(species, query, site, year, month, day, link) %>%
  summarise(num_weighed = length(weight))

me1 <- merge(bs_abund, unique(sp2[c("site", "year", "month", "day", "query", "link", "num_weighed")]), 
             by = c("site", "year", "month", "day", "query"), all.x = T)

grp1 <- filter(me1, count == num_weighed) #where #measured was number caught
grp2 <- filter(me1, count > num_weighed) #number measured is a subsample, need to expand (see below)
grp3 <- filter(me1, count < num_weighed) %>%#number measured is greater than number caught, problem here, probably go with measured as abundance
  mutate(count = num_weighed)
grp4 <- filter(me1, is.na(num_weighed)) #none measured, use entire dataset for weight distribution

grp1_3_s <- filter(sp1, link %in% grp1$link | link %in% grp3$link) %>%
  group_by(site, year, month, day, query, link) %>%
  summarise(mean = sum(weight)) %>%
  merge(me1[c("count", "link", "num_weighed")], ., by = "link", all.x = F, all.y = T) %>%
  mutate(count = num_weighed) %>%
  mutate(sd = NA)
  

#for each link draw from the weight distribution to the number in abundance;  resample distribution and average ####
#grp2 resample distribution
#grp2_s <- data.frame()

#for (i in unique(grp2$link)) {
#  t1 <- filter(grp2, link == i)
#  set.seed(12) 
#  t2 <- replicate(1000, {
#    t3 <- filter(sp1, link == i) %>%
#      pull(weight) %>%
#      sample(., t1$count, replace = T) %>%
#      sum()})
#  mean <- mean(t2)
#  sd <- sd(t2)
#  t4 = cbind(t1, mean, sd)
#  grp2_s = rbind(grp2_s, t4)
#}
#grp2_s <- grp2_s %>% select(-c("num"))

grp2_s <- sp1 %>% merge(., grp2[c("count", "link", "num_weighed")], by = "link", all.x = F, all.y = T) %>% 
  group_by(site, year, month, day, query, count,link, num_weighed) %>% 
  summarise(mean = mean(weight)) %>% 
  mutate(mean = mean * count,
         sd = NA)

#grp4_s draw from distributions for all surveys, take from gt4
distr_grp4 <- gt4 %>%
  filter(!is.na(weight)) # remove observations with missing weights

g1 <- grp4 %>%
  filter(query != "unla" & query != "unsc" & count != 0) 
g1$num <- as.character(as.numeric(factor(paste0(g1$species, g1$query, g1$site, g1$year, g1$month, g1$day))))
g1$link = paste("not_measured", g1$num, sep = "")

grp4_s <- data.frame()

for (i in unique(g1$link)) {
  p1 <- filter(g1, link == i)
  set.seed(12) 
  p2 <- replicate(100, {
    p3 <- filter(distr_grp4, query == p1$query) %>%
    pull(weight) %>%
    sample(., p1$count, replace = T) %>%
    sum()})
    mean <- mean(p2)
    sd <- sd(p2)
  p4 = cbind(p1, mean, sd)
  grp4_s = rbind(grp4_s, p4)
}
grp4_s <- grp4_s %>% select(-c("num"))

h1 <- rbind(grp1_3_s, grp2_s, grp4_s) %>%
  rename(weight = mean, abundance = count)


#make long data - fill in the gaps
h2 <- h1 %>%
  mutate(dat_site = as.character(factor(paste(site, year, month, day)))) %>%
  select(-c("link", "sd", "num_weighed")) %>%
  mutate(item_name = 1) %>%
  complete(dat_site, query, item_name) %>%
  separate(dat_site, c("site", "year", "month", "day"), sep = " ") %>%
  select(-c("item_name")) %>%
  replace(is.na(.), 0)


#calculate seine effort #####
seine_area <- read_csv("no_occupancy_model/Data/2022_10_31/bs_data/bs_seine_effort.csv") %>%
  filter(year >= 2018) %>%
  mutate(site = tolower(site))

seine_area <- seine_area %>%
 mutate(area = (3.1415 * (width1/2)^2) + (width1 * (extent1 - width1/2)) +   # set 1
                (3.1415 * (width2/2)^2) + (width2 * (extent2 - width2/2))) %>%  # set 2
   mutate(volume = (3.1415 * (width1/2)^2)*depth1 + (width1 * (extent1 - width1/2))*depth2/2 +   # set 1
                   (3.1415 * (width2/2)^2)*depth2 + (width2 * (extent2 - width2/2))*depth2/2)   # set 2
write_rds(seine_area, "no_occupancy_model/Data/2022_10_31/BS_effort.rds")


h3 <- filter(h2, year >= 2018)

eDNA_long <- read_rds("no_occupancy_model/Data/2022_10_31/eDNA_long.rds") #use to filter BS dataset

#calculate density metrics, clean up, filter
bs_long <- merge(h2[c("site", "year", "month", "day", "query", "abundance", "weight")], 
            seine_area[c("site", "year", "month", "day", "area", "volume")], 
            by = c("site", "year", "month", "day")) %>%
       filter(!is.na(weight))  %>%
  mutate(count_per_m2 = abundance/area)%>%
  mutate(count_per_m3 = abundance/volume)%>%
  mutate(grams_per_m2 = weight/area)%>%
  mutate(grams_per_m3 = weight/volume) %>%
  select(-c("volume", "area")) %>%
  mutate(dat_site = as.character(factor(paste(site, year, month, day)))) %>%
  mutate(date = make_date(year = year, month = month, day = day))%>%
  filter(dat_site %in% eDNA_long$dat_site) %>%
  mutate(p_a = if_else(abundance == 0, 0, 1))  #calculate presence absence

#remove taxa not detected in these surveys
u1 <- bs_long %>%
  select(c("query", "abundance")) %>%
  group_by(query) %>%
  summarise(total = sum(abundance)) %>%
  filter(total == 0)
bs_long <- bs_long %>% 
  filter(!query %in% u1$query)

write_rds(bs_long, "no_occupancy_model/Data/2022_10_31/bs_long.rds")







