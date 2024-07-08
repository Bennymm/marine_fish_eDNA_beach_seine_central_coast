#### Summarise "asv_by_site" from wrangling.R; count# of detections per site ####

#packages
library(tidyverse)
library(lubridate)

#data
dat <- readRDS("Data/2022_10_31/derived_data/pooled_asvMatrix.rds")
spec <- dat[c(2:1726)]                               #select Asv Matrix
out <- readRDS("Data/2022_10_31/derived_data/scratch/occProb_royallink.rds")
ASV_by_sample <- readRDS("Data/2022_10_31/derived_data/scratch/ASV_by_sample.rds")
rownames(out) <- out$X
out$X <- NULL
occprob <- out  

t1 <- ASV_by_sample %>%
  group_by(sample) 
t2 <- summarise_all(t1,list(sum))
t2 <- column_to_rownames(t2, "sample")
t3 <- t(t2)
t3 <- as.data.frame(t3)
sumreps <- t3

#### assign occupancy probability given # of detections ####
occprob_by_sample <- ifelse(sumreps == 3, occprob$PoO.3, ifelse(
  sumreps == 2, occprob$PoO.2, ifelse(
    sumreps == 1, occprob$PoO.1, occprob$PoO.0)))

#write.csv(occprob_by_sample, "Data/2022_10_31/derived_data/occprob_by_sample.csv")
saveRDS(occprob_by_sample, "Data/2022_10_31/derived_data/occprob_by_sample.rds")

# find all ASVs for which occprob is always less than 0.8 ... or 0.2
t9 <- ifelse(occprob_by_sample >= 0.8, 1, 0)
#write.csv(t9, "./Outputs/occupancy_modelling/royle_link/occ_by_sample.csv")
t10 <- as.data.frame(t(t9))
ASVcount_a <- as.data.frame(colSums(t10))           #find number of observations of each ASV
ASV0count_a <- filter(ASVcount_a, ASVcount_a[1] == 0)      #find ASVs with 0 observations
ASVswith0count_a <- rownames(ASV0count_a)
t10[,c(ASVswith0count_a)] <- NULL                   #remove ASVs with 0 "true positive" observations from ASV matrix
t15 <- as.matrix(t10)

#format ASV read abundance x site
t11 <- spec[c(colnames(t10))]                       #filter ASVs that have occupancy probability above threshold (defined above)
t12 <- cbind(dat$surveyID, t11)
colnames(t12)[1] <- "surveyID"                   #sum PCR read counts across samples
t13 <- t12 %>%
  group_by(surveyID) %>%
  summarise_all(list(sum)) %>%
  column_to_rownames("surveyID")
t16 <- as.matrix(t13)

ASV_LowOccRemoved <- as.data.frame(ifelse(t15 == 0, 0, t16)) # filter ASV read number matrix (t16) with occupancy thresolds (t15)
s1 <- ASV_LowOccRemoved %>%
  rownames_to_column(., var = "link") %>%
  pivot_longer(., !link, names_to = "ASV", values_to = "count")

site_survey_data <- read_rds("Data/2022_10_31/derived_data/sitesurvey_data.rds") %>%
  select(c("surveyID", "location", "date", "hakai_link_date"))

s2 <- merge(site_survey_data, s1, by.x = "surveyID", by.y = "link") %>%
  mutate(day = day(hakai_link_date), month = month(hakai_link_date), year = year(hakai_link_date)) %>%
  dplyr::rename(site = location) %>%
  mutate(dat_site = as.character(factor(paste(site, year, month, day)))) %>%
  select(c("ASV", "site", "surveyID", "dat_site",  "date", "hakai_link_date", "count")) %>%
  mutate(p_a = if_else(count == 0, 0, 1))                                                 #add presence/absence

#read index
s3 <- s2 %>%
 select(c("dat_site", "ASV", "count")) %>%
  pivot_wider(names_from = ASV, values_from = count) %>%
  column_to_rownames("dat_site")

#calculate read index  
PROP <- function(df) {sweep(df, MARGIN=2, colSums(df), FUN="/")}
#Take indices using a table with proportions of reads with families/OTUs in rows, samples in columns
INDEX <- function(df) {sweep(df, MARGIN=1, STATS=apply(df, MARGIN = 1, max), FUN="/")}

s4 <- s3 %>%
  t() %>% as.data.frame() %>%
  PROP() %>%
  INDEX() %>%
  t() %>% as.data.frame() 

eDNA_long <- s4 %>%
  rownames_to_column(var = "dat_site") %>%
  pivot_longer(!dat_site, names_to = "ASV", values_to = "index") %>%
  merge(s2,., by = c("dat_site", "ASV"))
  
saveRDS(eDNA_long, "Data/2022_10_31/derived_data/data12se_asvlong_lor_12s_ei.rds")
saveRDS(s4, "Data/2022_10_31/derived_data/data12se_asvmatrix_lor_12s_ei.rds")


#Taxonomy of removed ASVs

ASVswith0count_a

top10 <- read.delim("Data/2022_10_31/12S_ASV_sequences.length_var.blast.out.txt",
                    h=TRUE,
                    fill = TRUE)%>%
  `colnames<-`(c("ASV", "subject", "accesion_num", "taxa_ID", "perc_ID", "coverage", "evalue", "bitscore", "source", "taxonomy")) %>%
  as.data.frame() %>%
  na.exclude() %>%
  separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = " / ")

ASVs_removed_taxonomy <- filter(top10, ASV %in% ASVswith0count_a) %>%
  dplyr::select(c("ASV", "class")) %>%
  distinct()

summary_ASVs_removed_taxonomy <- ASVs_removed_taxonomy %>%
  group_by(class) %>%
  summarise(num = length(class))
776/842
842-776

observations_removed_ASV_by_sample <- ASV_by_sample[ASVswith0count_a]
number_of_hits <- as.data.frame(colSums(observations_removed_ASV_by_sample)) %>%
  `colnames<-`("hits")

summary_number_of_hits <- number_of_hits %>%
  group_by(hits) %>%
  summarise(num = length(hits))
737/842
842-737
