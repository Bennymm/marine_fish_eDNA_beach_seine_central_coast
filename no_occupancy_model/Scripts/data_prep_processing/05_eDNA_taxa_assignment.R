#read data ####
# NOTE: this will replace assignment_corrections04.R and eDNA_taxized05b_new.R 
# assign taxonomy
# create taxa_by_site_survey matrix
# plot some summaries

# packages and data ####

library(tidyverse)
library(RColorBrewer)
library(here)
library(vegan)
library(usedist)
library(taxize)
library(janitor)

ASVbysite <- readRDS("no_occupancy_model/Data/2022_10_31/derived_data/ASVmatrix_ei.rds")
sitesurvey_data <- readRDS("no_occupancy_model/Data/2022_10_31/derived_data/sitesurvey_data.rds")
sample_data <- read_csv("no_occupancy_model/Data/2022_10_31/Calvert_12S_metadata.csv")

top10 <- read.delim("no_occupancy_model/Data/2022_10_31/12S_ASV_sequences.length_var.blast.out20221101.txt",
                    h=TRUE,
                    fill = TRUE)

#format fish taxonomy
top10_fish <- top10 %>%
  `colnames<-`(c("ASV", "subject", "accesion_num", "taxa_ID", "perc_ID", "coverage", "evalue", "bitscore", "source", "taxonomy")) %>%
  as.data.frame() %>%
  na.exclude() %>%
  separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = " / ") %>%
  filter(class == "Actinopteri" | class == "Chondrichthyes") 

#select ASVs that were found in samples to be analyzed
top10_survey <- top10_fish %>%
  filter(ASV %in% colnames(ASVbysite))
        
#format taxonomy for non-fish
top10_notfish <- top10  %>%
  `colnames<-`(c("ASV", "subject", "accesion_num", "taxa_ID", "perc_ID", "coverage", "evalue", "bitscore", "source", "taxonomy")) %>%
  as.data.frame() %>%
  na.exclude() %>%
  separate(taxonomy, into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = " / ") %>%
  filter(class != "Actinopteri" & class != "Chondrichthyes")

#identify max percent ID for each ASV ####
max_ID <- top10_survey %>%
  group_by(ASV) %>%
  summarise(perc_ID = max(perc_ID))
#select taxonomy for max percent ID
f1 <- merge(max_ID, top10_survey, by = c("ASV", "perc_ID"))

#taxize #### only run once and turn off with "#" (as below) 
spec_unique <- unique(f1$species)
#check accepted naming and get higher taxonomy
top10_gbifid_higher <- classification(sci_id = spec_unique,
                                                   db = 'gbif',
                                                   #give back ID
                                                   return_id = TRUE) %>%
  #bind them together
  cbind(.) %>%
  rename("ids" = "query")
write_rds(top10_gbifid_higher, "no_occupancy_model/Data/2022_10_31/derived_data/top10_gbifid_higher.rds")
top10_gbifid_higher <- read_rds("no_occupancy_model/Data/2022_10_31/derived_data/top10_gbifid_higher.rds")

# merge with old names
new_taxonomy <- top10_gbifid_higher %>%
  mutate(old_species = spec_unique) %>%
  relocate(old_species, .after = species)
##################################manual editing required on next lines###########################################
#fix a few odd ducks manually, look over taxonomy assignment for missing species or odd names
new_taxonomy$species <- if_else(new_taxonomy$old_species == "Myzopsetta ferruginea", "Myzopsetta ferruginea", new_taxonomy$species) #missing name
new_taxonomy$genus <- if_else(new_taxonomy$old_species == "Myzopsetta ferruginea", "Myzopsetta", new_taxonomy$genus)    #missing name
new_taxonomy$species <- if_else(new_taxonomy$genus == "Mola", "Mola mola", new_taxonomy$species) #odd names


##################################manual editing required on next lines###########################################
write_csv(new_taxonomy, "no_occupancy_model/Data/2022_10_31/derived_data/taxonomyt2.csv")
#export and annotate: add "in_range" column in excell and annotate "y" if in 
#northeast pacific or tributaries, "n" if from other oceans (including western pacific)
#rename as below
new_taxa <- read_csv("no_occupancy_model/Data/2022_10_31/derived_data/taxonomy_range_annotated.csv") #import annotated taxa

#create dataframe with all top hits and accepted taxonomy ####
#select tazonomy for max percent ID
f2 <- merge(f1[c("ASV", "perc_ID", "species")], new_taxa[1:8], by = "species")

#count the numbers of families, genera, and species with equal max percent ID
f3 <- f2 %>%
  group_by(ASV, perc_ID) %>%
  summarise(across(c("family", "genus", "species"), ~ length(unique(.x)))) %>%
  `colnames<-`(c("ASV", "perc_ID", "fam_n", "gen_n", "spec_n"))
#group ASV to groups that need to be collapsed to family or within and genera or within
fam <- filter(f3, gen_n >= 2)
gen <- filter(f3, gen_n == 1 & spec_n > 1)
spec <- filter(f3, spec_n == 1)

#sort out multiple hits within family ####
#list groups where >2 species
fam_tax <- merge(fam, f2, by = c("ASV", "perc_ID")) %>%
  distinct()
#which groups have multiple species in range
r1 <- fam_tax %>%
  group_by(ASV) %>%
  summarize(n_in_range = sum(in_range == "y"))
r2 <- filter(r1, n_in_range == 1) 
d1 <- filter(fam_tax, ASV %in% r2$ASV) %>%      #ASVs with multiple hits, but only one species assignment in range, complete
  filter(in_range == "y") %>%
  mutate(LCT = species)%>%
  mutate(all_species = species) %>%
  mutate(level = "species") %>%
  .[c("ASV", "level", "LCT", "order", "family", "genus", "species", "all_species")]
r3 <- filter(r1, n_in_range > 1)                #when multiple species in range
r4 <- filter(fam_tax, ASV %in% r3$ASV) %>%
  filter(in_range == "y")
r5 <- r4 %>%
  group_by(ASV) %>%
  summarise(p1 = length(unique(genus))) %>%
  filter(p1 == 1)
r6 <- filter(r4, ASV %in% r5$ASV)               #ASVs with multiple species within genera, add to gen_tax below to cluster within genera
r7 <- filter(r4, !ASV %in% r5$ASV)              #ASVs with multiple genera within families, assign grouping manually
r8 <- r7[with(r7, order(ASV, species)), ] %>%
  group_by(ASV, order, family, in_range) %>%
  summarise(all_species = paste(species, collapse=", "))

##################################manual editing required on next lines###########################################
r9 <- data.frame(all_species = unique(r8$all_species))                         #table of family groups
r9$LCT <- c("Pleuronectidae1")    ######add a group name (in order) for each row in r9
r10 <- merge(r8, r9, by = "all_species") %>%
  mutate(level = "family") %>%
  mutate(genus = LCT) %>%
  mutate(species = LCT) %>%
  distinct()
d2 <- r10 %>%
  .[c("ASV", "level", "LCT", "order", "family", "genus", "species", "all_species")]

# merger q4, q5, q6 for table
tab_fam <- merge(r7[c("ASV", "species")], r8[c("ASV", "all_species")], by = "ASV") %>%
  merge(., r9, by = "all_species") %>%
  select(!ASV) %>%
  distinct()

#sort out multiple hits within genera ####
gen_tax <- merge(gen, f2, by = c("ASV", "perc_ID")) %>%
  distinct() %>%
  rbind(.,r6)                                   #add ASVs with multiple species within genera from above
#which groups have multiple species in range
q1 <- gen_tax %>%
  group_by(ASV) %>%
  summarize(n_in_range = sum(in_range == "y"))
q2 <- filter(q1, n_in_range == 1) 
d3 <- filter(gen_tax, ASV %in% q2$ASV) %>%      #ASVs with multiple hits, but only one species assignment in range, complete
  filter(in_range == "y")%>%
  mutate(LCT = species)%>%
  mutate(level = "species")%>%
  mutate(all_species = species) %>%
  .[c("ASV", "level", "LCT", "order", "family", "genus", "species", "all_species")]
q3 <- filter(q1, n_in_range > 1)                #when multiple species in range
q4 <- filter(gen_tax, ASV %in% q3$ASV) %>%
  filter(in_range == "y")                       #ASVs with multiple species in region, in a genus, assign to group below
q5 <- q4[with(q4, order(ASV, species)), ] %>%
  group_by(ASV, order, family, genus, in_range) %>%
  summarise(all_species = paste(species, collapse=", "))
q6 <- data.frame(all_species = unique(q5$all_species))                         #table of genus groups

##################################manual editing required on next lines########################################### 
q6$LCT <- c("Anoplarchus1", "Xiphister1", "Hexagrammos1", "Oligocottus1", "Pholis1", "Sebastes1", "Cottus1", ######add a group name (in order) for each row in q6#######
            "Myoxocephalus1", "Pholis1", "Gasterosteus aculeatus", "Lepidopsetta1")
q7 <- merge(q5, q6, by = "all_species") %>%
  mutate(level = "genus") %>%
  mutate(species = LCT)
d4 <- q7 %>%
  .[c("ASV", "level", "LCT", "order", "family", "genus", "species", "all_species")]
# merger q4, q5, q6 for table
tab_gen <- merge(q4[c("ASV", "species")], q5[c("ASV", "all_species")], by = "ASV") %>%
  merge(., q6, by = "all_species") %>%
  select(!ASV) %>%
  distinct()

##################################manual editing MAY BE required on next lines########################################### 
#species level issues: out of range ####
spec_tax <- merge(spec, f2, by = c("ASV", "perc_ID"))%>%
  distinct()
y1 <- spec_tax[c(1,6:13)] 
y2 <- y1 %>%                              #species outside of range, assign taxa by next best hit
  filter(in_range != "y") %>%
  mutate(LCT = c(""))%>%              ##### new assignment goes here (between "")#########
  mutate(level = c("")) %>%             ######### level of assignment goes here (between "") ###########
  mutate(in_range = c("y"))%>%         
  mutate(all_species = LCT)%>%         
  mutate(species = LCT)
y3 <- filter(y2, level == "genus") %>%
  mutate(all_species = species) %>%
  mutate(species = LCT)
y4 <- filter(y2, level == "family") %>%
  mutate(all_species = species) %>%
  mutate(species = species) %>%
  mutate(genus = LCT)
y5 <- filter(y2, level == "species") %>%
  mutate(all_species = species)%>%
  mutate(species = LCT) %>%
  mutate(genus = word(LCT,1))           #extract genus from species name
y6 <- filter(y2, level == "class") %>%
  mutate(all_species = species) %>%
  mutate(species = LCT) %>%
  mutate(genus = LCT) %>%
  mutate(family = LCT) %>%
  mutate(order = LCT) 
y7 <- rbind(y3,y4,y5,y6)      # NOTE some of y3-y6 don't do anything now,but may when we have other ID issues
d5 <- y7 %>%
  .[c("ASV", "level", "LCT", "order", "family", "genus", "species", "all_species")]


y8 <- y1 %>%                              #species inside of range, with single hits, or family or genus issues resolved by removing range issues
  filter(in_range == "y") %>%
  .[c("ASV","species",  "order", "family", "genus")] %>%
  rbind(.,d1[c("ASV","species",  "order", "family", "genus")]) %>%               #add species from family errors that were resolved by removing range issues
  rbind(.,d3[c("ASV","species",  "order", "family", "genus")])               #add species from genus errors that were resolved by removing range issues

# species that were assigned to groups for some ASVs but not others
y9 <- rbind(tab_fam, tab_gen) #grouping table
y10 <- merge(y8, y9, by = "species", all.x = T) 
y11 <- filter(y10, !is.na(y10$all_species)) %>%       #the ones that slipped through
  mutate(level = "genus")  %>%       #the ones that slipped through
  mutate(species = LCT)                                   #assign manually
y12 <- filter(y10, is.na(y10$all_species)) %>%
  mutate(level = "species") %>%
  mutate(LCT = species) %>%
  mutate(all_species = species)
d6 <- rbind(y12,y11)%>%
  .[c("ASV", "level", "LCT",  "order", "family", "genus", "species", "all_species")]

data <- rbind(d6, d5, d4, d2)

#summaries of ASVs ####

#####################################manual editing MAY BE required on next lines#########################
#There was a  problem with the one of the rockfish assignments (ASV4), is assigned to S. alutus, but a manual search
# beyond the top 10 shows 100s of other identical matches change here to Sebastes 2

data <- data %>%
  mutate(LCT = if_else(species == "Sebastes alutus", "Sebastes2", LCT)) %>%
  mutate(species = if_else(species == "Sebastes alutus", "Sebastes2", species)) %>%
  mutate(level = if_else(species == "Sebastes2", "genus", level)) %>%
  mutate(all_species = if_else(species == "Sebastes2", "Sebastes - alutianus, alutus, elongatus, flavidus, melanistictus, melanops,  ruberrimus", all_species)) %>%
  mutate(all_species = if_else(species == "Gasterosteus aculeatus", "Gasterosteus aculeatus", all_species)) %>%
  mutate(class = ifelse(species == "Hydrolagus colliei", "Chondrichthyes", "Actinopterygii"))
  


#number of ASVs with non-regional top hits
non_regional <- filter(f2, in_range == "n")
length(unique(non_regional$ASV)) # number of ASVs
unique(non_regional$genus) #genera

#occurences of taxa that were clustered to LITs
LIT_clusters <- filter(data, level != "species")
length(unique(LIT_clusters$ASV)) # number of ASVs
unique(LIT_clusters$genus) #genera
#27 species in LITs


#check that all ASVs are accounted for and not duplicated - they should all be the same
length(data$ASV)
length(unique(data$ASV)) 
length(unique(top10_survey$ASV))

table <- data[c("level", "LCT", "order", "family", "genus", "species", "all_species")] %>%
  distinct()
write_csv(table, "no_occupancy_model/Data/2022_10_31/derived_data/taxonomy_groups_12s_eDNA20221103.csv")

write_rds(data, "no_occupancy_model/Data/2022_10_31/derived_data/ASV_taxonomy_eDNA20221103.rds")


#make LCT by site matrix ####

LCT_by_site <- as.data.frame(t(ASVbysite)) %>%
  rownames_to_column( var = "ASV") %>%
  merge(data[c("ASV", "LCT")], ., by = "ASV", all.x = T) %>%
  select(-c("ASV")) %>%
  group_by(LCT) %>%
  summarise(across(.cols = everything(), ~sum(.x))) %>%
  column_to_rownames(var = "LCT") %>%
  t()%>%
  as.data.frame()

write_rds(LCT_by_site, "no_occupancy_model/Data/2022_10_31/LCT_by_site.rds")

#make table for augmenting Beach Seine taxonomy

b1 <- data
b2 <- b1 %>% separate(all_species, c("A", "B", "C", "D"), sep = ", ") ######### ignore warning ##########

b3 <- b2 %>% 
  select(-c("ASV")) %>%
  pivot_longer(cols = c("A", "B", "C", "D")) %>%
  rename('spec' = 'value') %>%
  select(-c("name")) %>%
  filter(!is.na(spec)) %>%
  distinct() %>%
  mutate(class = ifelse(species == "Hydrolagus colliei", "Chondrichthyes", "Actinopterygii"))

write_rds(b3, "no_occupancy_model/Data/2022_10_31/taxonomy_eDNA.rds")

#summaries of fish/non-fish detections ####

#select ASVs that passed occupancy models
top10_occ <- top10_fish %>%
  filter(ASV %in% colnames(ASVbysite))

#select ASVs that passed occupancy models
top10_occ_notfish <- top10_notfish %>%
  filter(ASV %in% colnames(ASVbysite))

length(unique(top10_survey$ASV)) #fish reads
length(unique(top10_notfish$ASV)) #non-target reads 
sum(LCT_by_site) #number of fish reads 
sum(ASVbysite) #total number of reads 
sum(ASVbysite) - sum(LCT_by_site)



sum(ASVbysite)


length(unique(top10_occ$ASV))
length(unique(top10_occ_notfish$ASV))

#select ASVs that passed occupancy models
top10_occ_notfish <- top10_notfish %>%
  filter(ASV %in% colnames(ASVbysite))

#identify max percent ID for each ASV ####
max_ID_notfish <- top10_occ_notfish %>%
  group_by(ASV) %>%
  summarise(perc_ID = max(perc_ID))
#select taxonomy for max percent ID
m1 <- merge(max_ID_notfish, top10_occ_notfish, by = c("ASV", "perc_ID"))





