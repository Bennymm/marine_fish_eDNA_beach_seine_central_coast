#Format environmental data, calculate temperature gradients and date difference, PCA of seawater mixing

#packages
library("vegan")
library("tidyverse")
library("dplyr")
library("here")
library("readxl")
library("sf")
library("rgdal")
library("lubridate")

#site survey dates
survey_dates <- read_rds("no_occupancy_model/Data/2022_10_31/eDNA_long.rds") %>%
  select(c("dat_site", "site", "date", "hakai_link_date")) %>%
  distinct() %>%
  mutate(dat_diff = abs(time_length(difftime(date, hakai_link_date), "days")))

#YSI data
ysi <- read_excel("no_occupancy_model/Data/2022_10_31/bs_data/FABSMasterData.xlsx", sheet = "ysi")
ysi_clean <- ysi %>%
  select(c("site", "year", "month", "day", "location", "temp")) %>%
  distinct() %>%                                                    #remove duplicate columns
  mutate(temp = if_else(temp == "NA", "", temp)) %>%
  pivot_wider(., names_from = location, values_from = temp)  %>%     #make wide
  select(!"NA") %>%
  mutate(a1 = as.numeric(a1)) %>%                                   #make numeric 
  mutate(b1 = as.numeric(b1)) %>%
  mutate(b2 = as.numeric(b2)) %>%
  mutate(c1 = as.numeric(c1)) %>%
  mutate(c2 = as.numeric(c2)) %>%
  mutate(c3 = as.numeric(c3)) %>%
  mutate(d1 = as.numeric(d1)) %>%
  mutate(d2 = as.numeric(d2)) %>%
  mutate(d3 = as.numeric(d3)) %>%
  mutate(d4 = as.numeric(d4)) %>%
  mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x)) %>%    #replace NA's with mean values
  mutate(temp.diff.b = abs(b1 - b2)) %>%                            #calculate temperature differences
  mutate(temp.diff.c = abs(c1 - c3)) %>%
  mutate(temp.diff.d = abs(d1 - d4)) %>%
  mutate(site = tolower(site)) %>%
  mutate(dat_site = paste(site, year, month, day, sep =  " "))

#write_rds(ysi_clean, "no_occupancy_model/Data/2022_10_31/bs_data/ysi.rds")
#ysi_clean <- readRDS("no_occupancy_model/Data/2022_10_31/bs_data/ysi.rds")

#habitat richness estimates (from rubidge) extracted in QGIS#
hab_rich <- readOGR("no_occupancy_model/Data/2022_10_31/bs_data/bs_2015_spatial_habrich.gpkg") #read geopackage
#format hab_rich
hab_rich <- as.data.frame(hab_rich)
hab_rich <- rename(hab_rich, site = field_1)
hab_rich <- hab_rich[c("site", "nFeat", "GiZScore", "GiPValue")]
write_rds(hab_rich, "no_occupancy_model/Data/2022_10_31/bs_data/hab_rich.rds")

#read in habitat and habitat distances
habitat <- read_csv("no_occupancy_model/Data/2022_10_31/bs_data/hakaiBS_habitat_modified.csv")
habitat $exposure_num <- as.integer(factor(habitat $exposure, 
                                      levels = c("vp", "p", "sp", "se", "e", "ve"), 
                                      labels = c(1,2,3,4,5,6)))

habitat$intertidal_primary_substrate_num <- as.integer(factor(habitat $intertidal_primary_substrate, 
                                                          levels = c("mud", "sand", "gravel", "cobble"), 
                                                          labels = c(1,2,3,4)))

habitat $subtidal_primary_substrate_num <- as.integer(factor(habitat $subtidal_primary_substrate, 
                                                        levels = c("mud", "sand", "gravel", "cobble"), 
                                                        labels = c(1,2,3,4)))

hab_dist1 <- read_csv("no_occupancy_model/Data/2022_10_31/bs_data/habitatdistance.csv")

h1 <- hab_dist1[c("site", "kelp", "seagrass", "water25m", "freshwater", "rockyshore")] %>%
  column_to_rownames("site")
h2 <- as.data.frame(ifelse(h1 >100, 0,1)) %>%
  mutate(h100m = rowSums(across(where(is.numeric)))) %>%
  select(h100m)
h3 <- as.data.frame(ifelse(h1 >200, 0,1)) %>%
  mutate(h200m = rowSums(across(where(is.numeric)))) %>%
  select(h200m)
h4 <- as.data.frame(ifelse(h1 >300, 0,1)) %>%
  mutate(h300m = rowSums(across(where(is.numeric)))) %>%
  select(h300m)
h5 <- as.data.frame(ifelse(h1 >400, 0,1)) %>%
  mutate(h400m = rowSums(across(where(is.numeric)))) %>%
  select(h400m)
h6 <- as.data.frame(ifelse(h1 >500, 0,1)) %>%
  mutate(h500m = rowSums(across(where(is.numeric)))) %>%
  select(h500m)
h7 <- as.data.frame(ifelse(h1 >1000, 0,1)) %>%
  mutate(h1000m = rowSums(across(where(is.numeric)))) %>%
  select(h1000m)
h8 <- as.data.frame(ifelse(h1 >2000, 0,1)) %>%
  mutate(h2000m = rowSums(across(where(is.numeric)))) %>%
  select(h2000m)

h9 <- cbind(h2, h3, h4, h5, h6, h7, h8) %>%
  rownames_to_column("site")
hab_dist <- merge( hab_dist1, h9, by = "site")  


#sediment data
sed <- read_csv("no_occupancy_model/Data/2022_10_31/bs_data/sediment.csv")

sed$silt <- sed$frac63um + sed$frac32um + sed$frac20um + sed$fracless20um
sed$silt_perc <- sed$silt/sed$sum
sed <- sed[c("site", "silt_perc")]
sed <- sed %>%
  group_by(site) %>%
  summarise(silt_percent = mean(silt_perc))

#merge datasets
# YSI - habitat - hab_rich - hab_dist - survey_dates - sed
#merge
m1 <- merge(ysi_clean, habitat, by = "site")
m2 <- merge(m1, hab_dist, by = "site")
m3 <- merge(m2, hab_rich, by = "site") %>%
  select(-c("site"))
m4 <- merge(survey_dates, m3, by = "dat_site")
m5 <- merge(m4,sed, by = "site")

write_rds(m5, "no_occupancy_model/Data/2022_10_31/environmental.rds")

