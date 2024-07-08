### Goal = convert read numbers to eDNA index

## read in packages
library(tidyverse)
library(vegan)
library(here)

ASVbysite <- readRDS("no_occupancy_model/Data/2022_10_31/derived_data/ASVmatrix.rds") %>%
  column_to_rownames("surveyID") 
ASVbysite = ASVbysite[,colSums(ASVbysite) > 0]
colSums(ASVbysite)  

##read in functions
#index from (Jacobs-Palmer et. al., 2020)
#Take proportions from table of raw read counts with families/OTUs in rows, samples in columns
PROP <- function(df) {sweep(df, MARGIN=2, colSums(df), FUN="/")}
#Take indices using a table with proportions of reads with families/OTUs in rows, samples in columns
INDEX <- function(df) {sweep(df, MARGIN=1, STATS=apply(df, MARGIN = 1, max), FUN="/")}

t1 <- PROP(ASVbysite)
t2 <- INDEX(t1)

saveRDS(t2, "no_occupancy_model/Data/2022_10_31/derived_data/ASVmatrix_ei.rds")


