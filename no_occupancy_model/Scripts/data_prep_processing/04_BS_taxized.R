#assign correct(ed) taxonomy to fish captured in beach seines

# packages ####
library(here)
library(taxize)
library(tidyverse)
library(janitor)

taxaBS <- read.csv("no_occupancy_model/Data/2022_10_31/bs_data/fishcodes.csv")

colnames(taxaBS)[3] <- "species"
colnames(taxaBS)[1] <- "query"
taxaBS <- filter(taxaBS, common_name != "no catch")

## Beach seine taxonomy clean up data content ####
taxaBS <- taxaBS %>%
  # remove anything that contains "unknown"
  mutate(across(.cols = everything(),
                .fns = ~ str_replace_all(string = .,
                                         pattern = "unknown.*",
                                         replacement = NA_character_))
  ) %>%
  # remove "no identification"
  mutate(across(.cols = everything(),
                .fns = ~ str_replace_all(string = .,
                                         pattern = "^no identification$",
                                         replacement = NA_character_))
  ) %>%
  # remove anything that contains "uncultured"
  mutate(across(.cols = everything(),
                .fns = ~ str_replace_all(string = .,
                                         pattern = "uncultured.*",
                                         replacement = NA_character_))
  ) %>%
  # filter out "filtered out" hits
  filter(species != "filtered out")

## clean up taxonomic assignment ----
#strip higher taxonomy and repopulate with correct assignments
# resolve species names
data12se_species_resolved <- tol_resolve(taxaBS %>%
                                           #remove nas
                                           drop_na(species) %>%
                                           #take species vector
                                           pull(species) %>% #length() - 606 asvs identified to species
                                           #take distinct values
                                           unique(.) #%>% length() - 87 unique species identified
)

data12se_species_resolved %>%
  pull(is_synonym) %>%
  unique() #none are synonyms so use original names

# get gbif ids for the scientific names
data12se_species_gbifid <- get_gbifid(sci = taxaBS %>%
                                        #remove nas
                                        drop_na(species) %>%
                                        #take species vector
                                        pull(species) %>% #length() 
                                        #take distinct values
                                        unique(.), #%>% length()
                                      #run interactively - if more than one id is found, ask me to chose - check taxonomy_decisions.txt for decisions made
                                      # input 1 each time it asks
                                      ask = T) #found 153/153 gbifids
#make into a dataframe
data12se_species_gbifid_df <- data12se_species_gbifid %>% 
  as.data.frame(.) %>%
  #add in species info
  mutate(species = taxaBS %>%
           #remove nas
           drop_na(species) %>%
           #take species vector
           pull(species) %>% #length() 
           #take distinct values
           unique(.)
  )

#get higher taxonomy ####
data12se_species_gbifid_higher <- classification(sci_id = data12se_species_gbifid_df %>%
                                                   #filter out the ones not found
                                                   filter(match == "found") %>%
                                                   pull(ids) %>%
                                                   unique(), #%>% length(), 
                                                 db = 'gbif',
                                                 #give back ID
                                                 return_id = TRUE) %>%
  #bind them together
  cbind(.) %>%
  mutate(class = "Actinopterygii")


data12se_species_gbifid_higher %>% pull(species) %>% unique() %>% length() 

# reformat table
data12se_species_gbifid_higher <- data12se_species_gbifid_higher %>% 
  # rename query row
  rename("gbif_query" = "query") %>% 
  #add in asv associated with each species - right join because there can be more than one
  right_join(., 
             taxaBS %>%
               #take columns of interest
               select(species,
                      query) %>%
               #remove rows containing NA
               drop_na() %>%
               #merge with gbifid
               full_join(.,
                         data12se_species_gbifid_df %>% 
                           #select columns of interest
                           select(species,
                                  ids) %>%
                           #rename column for later merge
                           rename("gbif_query" = "ids")
               ) %>%
               #take columns for merge
               select(gbif_query, 
                      query),
             by = "gbif_query") %>%
  #drop rows of species not found by gbif
 # drop_na(species_id) %>%
  #add indicator for level of query
  mutate(taxo_query_level = "species") %>%
  #add query source database
  mutate(higher_taxo_source = "gbif") %>%
  #clean column names
  clean_names()

data12se_species_gbifid_higher %>% pull(species) %>% unique() %>% length() 
data12se_species_gbifid_higher %>% pull(query) %>% unique() %>% length() 

saveRDS(data12se_species_gbifid_higher,"no_occupancy_model/Data/2022_10_31/derived_data/taxa_BS_resolved.rds")


