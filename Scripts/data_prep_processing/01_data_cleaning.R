### Initial Cleaning of eDNA data ###

#Packages ####
library(tidyverse)
library(here)
library(vegan)
library(usedist)
library(lubridate)
library(GUniFrac)

#read data ####
asv_matrix <- read.table("Data/2022_10_31/sequence_table.12S.merged.w_ASV_names.length_var.merged_dataset.txt",
                         header = T)

sample_data <- read_csv("Data/2022_10_31/Calvert_12S_metadata.csv")
taxa <- read.delim("Data/2022_10_31/12S_ASV_sequences.length_var.blast.out.txt",
                   h=TRUE,
                   fill = TRUE)

#formatting and partitioning ####
#rename columns
colnames(sample_data) <- c("SampleID_rep", "sample_ID", "PCR_rep", "sample_type", "location", "hakai_project", "date",
                           "hakai_link_date", "year", "depth", "spatial_rep", "location_name_old", "location_name",
                           "experiment", "temp_surf", "sal_surf", "pH_surf",
                           "temp_bot", "sal_bot", "pH_bot", "bottom_depth", "vol_filtered",
                           "lab_notes")
colnames(asv_matrix)[1] <- "SampleID_rep"
#fix experiment column
sample_data$experiment <- str_remove(sample_data$experiment, "repeat-")
sample_data$experiment <- ifelse(sample_data$experiment == 0, "surface", sample_data$experiment)
#merge sample data and matrix
sampledata_asvmatrix <- merge(sample_data, asv_matrix, by = "SampleID_rep")
#remove blanks, make DF of blanks
controls <- filter(sampledata_asvmatrix, sample_type != "biological")
sampledata_asvmatrix <- filter(sampledata_asvmatrix, sample_type == "biological" & experiment != "grid")

b1 <- gather(sampledata_asvmatrix, ASV, count, c(24:1748))
b2 <- gather(controls, ASV, count, c(24:1748))

b3<- b2 %>%
  group_by(sample_type, ASV) %>%
  summarise(mean_count = mean(count))

ext_blank <- filter(b3, sample_type == "Extraction_blank")
PCRneg_2018 <- filter(b3, sample_type == "PCR_Negative_2018")
PCRneg_2019_2020 <- filter(b3, sample_type == "PCR_Negative_2019_2020")

b4 <- merge(b1, ext_blank[c("ASV", "mean_count")], by=c("ASV")) %>% rename("ext_bl_count" = "mean_count")

b2018 <- filter(b4, year == 2018) %>%
  merge(., PCRneg_2018[c("ASV", "mean_count")], by=c("ASV")) %>% rename("PCRneg_count" = "mean_count")

b20192020 <- filter(b4, year == 2019 | year == 2020)

#pool PCR reps
t1 <- sampledata_asvmatrix[c(2,24:1748)]
pooled_asvmatrix <- t1 %>%
  group_by(sample_ID) %>%
  summarise_all(~sum(.x))
t2 <- sampledata_asvmatrix[c(2,4:23)]
t3 <- as.data.frame(t2[!duplicated(t2), ])
t4 <- merge(t3, pooled_asvmatrix, by = "sample_ID")
t4$loc_dat_depth <- as.numeric(factor(paste0(t4$location, t4$date, t4$experiment)))

t4 <- t4 %>% relocate(loc_dat_depth, .after = location)
t5 <- t4 %>%
  mutate(surveyID = paste("ss", loc_dat_depth, sep = ""))%>%
  relocate(surveyID, .after = location)

pooled_asvmatrix <- t5[c(4,24:1748)] #select columns for sample ID and all ASVs
sample_data <- t5[1:23]

#% similarity PCR reps  ####
#use pooled_asvmatrix grouping on site_survey
#calculate bray curtis dissimilarity
a1 <- vegdist(pooled_asvmatrix[2:764], method = "bray", binary = T)

# calculate distance from group centroid using Multivariate homogeneity of groups dispersions (variances)
a2 <- betadisper(a1, group = pooled_asvmatrix$surveyID, type = "centroid")

#plot with 95% CI elipses
plot(a2, hull = F, ellipse = T, conf = 0.95, label = F, main = "95% CI elipse around sample replicates")

#extract distances
a3 <- cbind(pooled_asvmatrix$surveyID, as.data.frame(a2$distances))
colnames(a3) <- c("sample_name", "distances")

#calculate 95% CIs
a4 <- a3 %>%
  group_by(sample_name) %>%
  mutate(meandis=mean(distances), # mutate when grouped applies the summary to each row
         maxdist=max(distances),
         SDdis=sd(distances),
         CIdis = 1.96*SDdis/sqrt(3), 
         upperCIdis=meandis+CIdis) #calculate upper end of 95% CI

max(a4$meandis) #no dissimilarity is greater than 0.49 (Kelly 2018)

#check that PCR reps are within CIs
a4$diff <- ifelse(a4$upperCIdis<a4$distances, "problem", "ok")
a5 <- filter(a4, diff == "problem")

#collapse sample data to site ####

site_data <- sample_data
  #summarize temp, depth, pH and volume
env_mean <- sample_data %>%
  group_by(surveyID) %>%
  summarise_at(vars(temp_surf, sal_surf ,pH_surf, temp_bot, sal_bot, pH_bot), .funs = mean, na.rm = TRUE) 
depth_max <- sample_data %>% #because some data entered as "0" for stations where YSI wasn't done
  group_by(surveyID) %>%
  summarise_at(vars(bottom_depth), .funs = max, na.rm = TRUE) 
vol_total <- sample_data %>% #because some data entered as "0" for stations where YSI wasn't done
  group_by(surveyID) %>%
  summarise_at(vars(vol_filtered), .funs = sum, na.rm = TRUE) 
#remove above vars
site_data[c("temp_surf", "sal_surf", "pH_surf", "temp_bot", "sal_bot", "pH_bot", 
              "bottom_depth", "vol_filtered", "lab_notes", "sample_ID", "spatial_rep",
              "location_name_old", "location_name", "depth")] <- NULL
  #remove duplicates
site_data <- as.data.frame(site_data[!duplicated(site_data), ])
  #merge back together
site_data <- site_data %>% 
  merge(.,env_mean, by = "surveyID") %>% 
  merge(.,depth_max, by = "surveyID") %>% 
  merge(.,vol_total, by = "surveyID") 

site_data$location <- tolower(site_data$location)

#Rarefaction ####
#identify minimum sequence depth - many sites can be discarded that aren't used in analyses
u1 <- data.frame(read_depth = rowSums(pooled_asvmatrix[2:ncol(pooled_asvmatrix)]), sample = pooled_asvmatrix$surveyID)
#23988 is the minimum read depth for a paired eDNA-beach seine
u2 <- filter(u1, read_depth < 23987) %>% # only one is paired with beach seine
  select(c("sample")) %>%
  distinct()

#select site-surveys to retain
k1 <- sample_data %>%
  filter(experiment == "surface") %>%         #surface samples
  filter(hakai_project == "beach_seine") %>%  #BS
  filter(!is.na(hakai_link_date)) %>%         #has a paired survey
  filter( !surveyID %in% u2$sample)
k2 <- filter(pooled_asvmatrix, surveyID %in% k1$surveyID)
sum(k2[c(2:ncol(k2))])  #24629880 reads across all sites, before rarefying

#rarefy or proportional
o2 <- Rarefy(k2[2:ncol(k2)], 23988)
o3 <- as.data.frame(o2$otu.tab.rff) %>%
  mutate(surveyID = k2$surveyID) %>%
  relocate(surveyID)
sum(o3[c(2:ncol(o3))])  #3310344 read after rarefying

#before rarefication summaries
R1 <- k2 %>%
 dplyr::select(-c("surveyID")) %>%
  colSums() %>%
  as.data.frame() %>%
  `colnames<-`("reads") %>%
  filter(reads > 0)
R1_01 <- ifelse(k2 == 0, 0,1) %>%
  as.data.frame() %>%
  dplyr::select(-c("surveyID")) 
length(R1$reads) #number of ASVs
sum(R1$reads)  #number of reads
sum(R1_01)  #number of observations
#after rarefication summaries
R2 <- o3 %>%
  dplyr::select(-c("surveyID")) %>%
  colSums() %>%
  as.data.frame() %>%
  `colnames<-`("reads") %>%
  filter(reads > 0)
R2_01 <- ifelse(o3 == 0, 0,1) %>%
  as.data.frame() %>%
  dplyr::select(-c("surveyID")) 
length(R2$reads) #number of ASVs
sum(R2$reads)   #number of reads
sum(R2_01)  #number of observations

4690-4638
rownames(R1)[!rownames(R1) %in% rownames(R2)] # removed by rarefaction


#save new RDS files ####
saveRDS(site_data, file = "Data/2022_10_31/derived_data/sitesurvey_data.rds")
saveRDS(sample_data, file = "Data/2022_10_31/derived_data/sample_data.rds")
saveRDS(o3, file = "Data/2022_10_31/derived_data/pooled_asvMatrix.rds")
saveRDS(pooled_asvmatrix, file = "Data/2022_10_31/derived_data/pooled_asvMatrix_no_rar.rds")
