#Prepare and run occupancy models for each ASV

#packages
library(tidyverse)
#library(jagsUI)
library(here)

#data
dat <- readRDS("Data/2022_10_31/derived_data/pooled_asvMatrix.rds")
spec <- dat[c(2:1726)]                               #select ASV Matrix
spec01 <- as.data.frame(ifelse(spec == 0, 0, 1))     #convert to binary

#summaries and cleaning
ASVcount <- as.data.frame(colSums(spec01))           #find number of observations of each ASV
ASVcountNo0s <- filter(ASVcount, ASVcount[1] != 0)   #count of detections by ASV
ASV0count <- filter(ASVcount, ASVcount[1] == 0)      #find ASVs with 0 observations
ASVswith0count <- rownames(ASV0count)
ASV1count <- filter(ASVcount, ASVcount[1] == 1)      #find ASVs with 1 observation (one PCR detection in dataset)
ASVswith1count <- rownames(ASV1count)
spec[,c(ASVswith0count)] <- NULL                     #remove ASVs with 0 observations from ASV matrix
spec01[,c(ASVswith0count)] <- NULL                   #remove ASVs with 0 observations from ASV matrix
ASV_by_sample <- cbind(dat$surveyID, spec01)      #join with sample name
colnames(ASV_by_sample)[1] <- "sample"               #rename column

# ASVs and read counts if just removing singletons (a hypothetical; not to be used further)
spec_s <- spec
spec_s01 <- spec01
spec_s[,c(ASVswith1count)] <- NULL                     #remove ASVs with 0 observations from ASV matrix
spec_s01[,c(ASVswith1count)] <- NULL 
sum(spec_s)

saveRDS(ASV_by_sample, "Data/2022_10_31/derived_data/scratch/ASV_by_sample.rds")

#formatt for occupancy models
#make a list of dataframes, if you get "Error: `n()` must only be used inside dplyr verbs." restart R. 
# There is a conflict with one of the packages from OccupancyModel.R

ASVs <- colnames(ASV_by_sample)
ASVs <- ASVs[-1]
ASVlist <- list()

for(i in ASVs){
  t1 <- dplyr::select(ASV_by_sample, sample, i)
  ASV <- t1 %>% 
    group_by(sample) %>%
    mutate(id = paste0("X", 1:n())) %>%
    spread(id, i)
  ASV$sample <- NULL
  ASVlist[[length(ASVlist)+1]] = ASV
}

ASVlist[7] #check to see that it is formatted properly - should be 3 columns X1,X2,X3 and 1/0 in rows

saveRDS(ASVlist, file="Data/2022_10_31/derived_data/scratch/ASVlist.rds")
saveRDS(ASVs, file="Data/2022_10_31/derived_data/scratch/ASVs.rds")


