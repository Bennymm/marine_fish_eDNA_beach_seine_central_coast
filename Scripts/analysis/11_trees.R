#Build a dendrogram from higher taxonomy and colour code by method. 
#Plot frequency of positive detections by method.

# packages and loading data ####
library(ape)
library(tidyverse)
#BiocManager::install("ggtree")
library(ggtree)
library(phylobase)

#Data
long <- read_rds("Data/2022_10_31/derived_data/shared_long.rds")

taxab <- read_rds("Data/2022_10_31/derived_data/BS_taxonomy_reconciled.rds")%>%
  mutate(genus = if_else(LCT_shared == "Pleuronectidae1", "Pleuronectidae1", genus))%>%
  mutate(genus = if_else(genus == "Pleuronectidae1", "Pleuronectidae_", genus)) # a few stragglers from formatting these tables

taxae <- read_rds("Data/2022_10_31/derived_data/ASV_taxonomy_reconciled.rds")%>% 
  mutate(LCT_shared = if_else(LCT_shared == "Sebastes alutus", "Sebastes3", LCT_shared))%>%
  mutate(genus = if_else(genus == "Pleuronectidae1", "Pleuronectidae_", genus)) # a few stragglers from formatting these tables

all <- rbind(taxab[c("class", "order", "family", "genus", "LCT_shared")],
             taxae[c("class", "order", "family", "genus", "LCT_shared")]) %>%
  distinct() %>%
  mutate(genus = if_else(genus == "Pleuronectidae1", "Pleuronectidae_", genus)) # a few stragglers from formatting these tables



#determine shared/exclusive detections across taxonomic levels
ssh <- taxae$LCT_shared[taxae$LCT_shared %in% taxab$LCT_shared]
se <- taxae$LCT_shared[!taxae$LCT_shared %in% taxab$LCT_shared]
sb <- taxab$LCT_shared[!taxab$LCT_shared %in% taxae$LCT_shared]

gsh <- taxae$genus[taxae$genus %in% taxab$genus]
ge <- taxae$genus[!taxae$genus %in% taxab$genus]
gb <- taxab$genus[!taxab$genus %in% taxae$genus]

fsh <- taxae$family[taxae$family %in% taxab$family]
fe <- taxae$family[!taxae$family %in% taxab$family]
fb <- taxab$family[!taxab$family %in% taxae$family]

osh <- taxae$order[taxae$order %in% taxab$order]
oe <- taxae$order[!taxae$order %in% taxab$order]
ob <- taxab$order[!taxab$order %in% taxae$order]

csh <- taxae$class[taxae$class %in% taxab$class]
ce <- taxae$class[!taxae$class %in% taxab$class]
cb <- taxab$class[!taxab$class %in% taxae$class]

psh <- taxae$phylum[taxae$phylum %in% taxab$phylum]
pe <- taxae$phylum[!taxae$phylum %in% taxab$phylum]
pb <- taxab$phylum[!taxab$phylum %in% taxae$phylum]

sh <- c(ssh, gsh, fsh, osh, csh, psh) %>%
  na.omit(.) %>%
  .[!duplicated(.)]%>%
  as.data.frame() %>%
  mutate(group = "detected with both methods")

e <- c(se, ge, fe, oe, ce, pe) %>%
  na.omit(.) %>%
  .[!duplicated(.)] %>%
  as.data.frame() %>%
  mutate(group = "detected with eDNA only")

b <- c(sb, gb, fb, ob, cb, pb) %>%
  na.omit(.) %>%
  .[!duplicated(.)]%>%
  as.data.frame() %>%
  mutate(group = "detected with beach seine only")

taxa_groups <- rbind(sh, e, b) %>%
  rename(., taxa = .) %>%
  distinct()

s1 <- as.data.frame(ssh) %>%
  rename(taxa = ssh) %>%
  mutate(group = "detected with both methods")
e1 <- as.data.frame(se) %>%
  rename(taxa = se) %>%
  mutate(group = "detected with eDNA only")
b1 <- as.data.frame(sb) %>%
  rename(taxa = sb) %>%
  mutate(group = "detected with beach seine only")

spec_shared <- rbind(s1,e1,b1) %>%
  distinct()

#dendrogram with higher taxonomy ####

t1 <- as.data.frame(lapply(all, factor))
t2 <- as.phylo(~class/order/family/genus/LCT_shared, data = t1, collapse=FALSE)
t3 <- as(t2, 'phylo4')
t4 <- as.data.frame(t3@label) %>%
  rename(taxa = "t3@label") %>%
  merge(., spec_shared, by = "taxa") %>%
  column_to_rownames("taxa")
t5 <- phylo4d(t3, t4)

t6 <- as.data.frame(t3@label) %>%
  rename(taxa = "t3@label") %>%                         
  rownames_to_column() %>%
  merge(., taxa_groups, by = "taxa", all.x = T) %>%
  column_to_rownames("rowname")

t6 <- as.data.frame(t3@label) %>%
  rename(taxa = "t3@label") %>%
  merge(., taxa_groups, by = "taxa") %>%
  column_to_rownames("taxa")

#make a tree and annotate
tree_v1 <-  ggtree(t5, open.angle = T, size = 0.5, aes(color=group)) +
  scale_color_manual(name = "Observation type",
                     values=c("detected with beach seine only" = "#E69F00",
                              "detected with eDNA only" = "#D55E00",
                              "detected with both methods" = "#000000")) +
  geom_tiplab(size = 3, as_ylab=F, offset = 0.2) +
  theme(legend.position = c(0.24,0.95),
        legend.background = element_rect(fill = "transparent"),
        legend.key = element_rect(colour = NA, fill = NA),
        plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
        aspect.ratio = 2) +
  coord_cartesian(xlim = c(0, 6), 
                          clip = 'off') +
  annotate("text", x = 2.5, y= 73.9, label = "Sculpins", angle = 0, size = 2.6) +
  annotate("text", x = 3, y= 64.6, label = "Rockfishes", angle = 0, size = 2.6) +
  annotate("text", x = 2.5, y= 58.9, label = "Greenlings", angle = 0, size = 2.6) +
  annotate("text", x = 2.5, y= 54.6, label = "Snailfishes", angle = 0, size = 2.6) +
  annotate("text", x = 2.5, y= 52.1, label = "Poachers", angle = 0, size = 2.6) +
  annotate("text", x = 2.5, y= 50.1, label = "Hemitripteridae", angle = 0, size = 2.6) +
  annotate("text", x = 2.5, y= 46.5, label = "Perches", angle = 0, size = 2.6) +
  annotate("text", x = 2.5, y= 42.1, label = "Pricklebacks", angle = 0, size = 2.6) +
  annotate("text", x = 2.5, y= 38.6, label = "Gobies", angle = 0, size = 2.6) +
  annotate("text", x = 3, y= 36.1, label = "Kelpfishes", angle = 0, size = 2.6) +
  annotate("text", x = 2.5, y= 34.1, label = "Gunnels", angle = 0, size = 2.6) +
  annotate("text", x = 1.5, y= 17.4, label = "Flatfishes", angle = 0, size = 2.6) +
  annotate("text", x = 2, y= 24.4, label = "Salmonids", angle = 0, size = 2.6) +
  annotate("text", x = 1.5, y= 12.5, label = "Cods", angle = 0, size = 2.6) +
  annotate("text", x = 2, y= 10.1, label = "Clingfishes", angle = 0, size = 2.6) +
  annotate("text", x = 1.5, y= 8.1, label = "Clupeiformes", angle = 0, size = 2.6) +
  annotate("text", x = 1.5, y= 6.1, label = "Gasterosteiformes", angle = 0, size = 2.6)

# geom_text(aes(label=node), hjust = -0.2) #use to diagnose
tree_v1

#flip a few branches to aid in discussion
#nt <- flip(tree_v1, 48, 47)

ggsave("./Figures_tables/tree_by_method.png", 
       plot = tree_v1,
       width = 20, height = 32, units = "cm")



#create stacked bar chart of prevalence for each LCT in tree (in order)

#get tip labels in order 
extract_tree_data <- function(tree_disp, displayorder=TRUE) {
  td_out <- tree_v1$data
  if (displayorder) {
    td_out <- dplyr::arrange(td_out,y)
  }
  return(td_out)
}
LCTs_in_order <- extract_tree_data(ggt_reordered) %>% 
  dplyr::filter(isTip) %>% 
  dplyr::pull(label) 
write_rds(LCTs_in_order, "Data/2022_10_31/derived_data/LCTs_in_order.rds")

#make counts of detection error type

q1 <- long %>%
  select(c("LCT_shared", "p_a_bs", "p_a_eDNA")) %>%
  mutate(shared_p = if_else(p_a_bs == 1 & p_a_eDNA == 1, 1, 0)) %>%
  mutate(eDNA_p = if_else(p_a_bs == 0 & p_a_eDNA == 1, 1, 0)) %>%
  mutate(bs_p = if_else(p_a_bs == 1 & p_a_eDNA == 0, 1, 0)) %>%
  mutate(shared_a = if_else(p_a_bs == 0 & p_a_eDNA == 0, 1, 0)) %>%
  group_by(LCT_shared) %>%
  summarise(across(c("shared_p", "eDNA_p" , "bs_p", "shared_a"), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(!LCT_shared, names_to = "error", values_to = "count")

q1$group1 <- factor(q1$error, levels = c("bs_p",
                                         "shared_a", #assign factors
                                         "shared_p",
                                         "eDNA_p"))
q1 <- with(q1, q1[order(group1),])   

#plot
frequency <- 
ggplot(q1, aes(y=count, x=LCT_shared, fill = group1)) + 
  geom_bar(position="stack", stat="identity", width = 0.5) +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust = 0.2, size = 8.5)) +
  scale_x_discrete(limits = LCTs_in_order) +
  theme_classic() +
#  scale_color_manual(values=c("#000000", "#E69F00", "#D55E00", "light grey")) +
  scale_fill_manual(values=c("#E69F00", "light grey", "#000000","#D55E00")) +
  ylab("Number of observation types")
  
frequency 
ggsave("./Figures_tables/detection_frequency.png", 
       plot = frequency,
       width = 13, height = 32, units = "cm")

