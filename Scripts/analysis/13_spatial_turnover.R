
#plot spatial pairwise dissimilarities for each method
#model turnover as a function of environmental variables and pairwise distance
# packages ####
library(tidyverse)
library(here)
library(lme4)
library(vegan)
library(ape)
library(RColorBrewer)
library(gdm) #https://cran.r-project.org/web/packages/gdm/gdm.pdf
library(dplyr)
library(ade4)
library(ggpubr)
library(GUniFrac)

#data and prep ####
#long data to calculate richness
long <- read_rds("Data/2022_10_31/derived_data/shared_long.rds")
#long <- filter(long, gamma_error == "both")             # run this to do analysis for taxa detected in both methods
#envonmental data
env <- read_rds("Data/2022_10_31/derived_data/environmental.rds") %>%
  mutate(exposure_num = ifelse(exposure_num == 5,4,exposure_num))

eDNA_mat <- pivot_wider(long[c("LCT_shared", "dat_site", "p_a_eDNA")], 
                  names_from = LCT_shared,
                  values_from = p_a_eDNA) %>%
  column_to_rownames(var = "dat_site")
eDNA_mat <- eDNA_mat[, which(colSums(eDNA_mat) != 0)]
eDNA_mat <- eDNA_mat[ order(row.names(eDNA_mat)), ]

BS_mat <- pivot_wider(long[c("LCT_shared", "dat_site", "p_a_bs")], 
                  names_from = LCT_shared,
                  values_from = p_a_bs) %>%
  column_to_rownames(var = "dat_site")
BS_mat <- BS_mat[, which(colSums(BS_mat) != 0)]
BS_mat <- BS_mat[ order(row.names(BS_mat)), ]

BS_mat_count <- pivot_wider(long[c("LCT_shared", "dat_site", "p_a_bs")], 
                      names_from = LCT_shared,
                      values_from = p_a_bs) %>%
  column_to_rownames(var = "dat_site")
BS_mat_count <- BS_mat_count[, which(colSums(BS_mat_count) != 0)]
BS_mat_count <- BS_mat_count[ order(row.names(BS_mat_count)), ]

o2 <- rrarefy(BS_mat, 700) 

shared_species <- colnames(eDNA_mat)[colnames(eDNA_mat) %in% colnames(BS_mat)]     #turn on to run analysis with only shared species
eDNA_mat <- eDNA_mat %>% dplyr::select(c(shared_species))
BS_mat <- BS_mat %>% dplyr::select(c(shared_species))

#Plot components of dissimilarity by method ####

#turnover
turn_e <- designdist(eDNA_mat, "2 * pmin(b,c) / (a + 2 * pmin(b,c))", abcd = T)
turn_b <- designdist(BS_mat, "2 * pmin(b,c) / (a + 2 * pmin(b,c))", abcd = T)

p1 <- as.data.frame(matrix(turn_e, dimnames=list(t(outer(colnames(turn_e), rownames(turn_e), FUN=paste)), NULL)))
p2 <- as.data.frame(matrix(turn_b, dimnames=list(t(outer(colnames(turn_b), rownames(turn_b), FUN=paste)), NULL)))

dat_turn <- data.frame(turn_e = p1, turn_b = p2) %>% # a better way to make a dataframe ! don't use tibble()
  `colnames<-` (c("turn_e","turn_b"))

plot_turn <- ggplot(data = dat_turn, aes(x = turn_b, y = turn_e)) +
  geom_point(size = 0.3) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1)) + 
  scale_y_continuous(limits = c(0, 1)) +
  annotate("segment", x = 0, xend = 0.9, y = 0, yend = 0.9) +
  annotate("text", x = 0.95, y = 0.95, label = "1:1") +
  ylab("Pairwise turnover - eDNA") +
  xlab("Pairwise turnover - beach seining") +
  annotate("text", x = 0.01, y = 1, label = "A")
plot_turn
ggsave("./Figures_tables/spatial_turnover.png", 
       plot = plot_turn,
       width = 3, height = 3, units = "in")

mean(dat_turn$turn_e)
sd(dat_turn$turn_e)
mean(dat_turn$turn_b)
sd(dat_turn$turn_b)
t.test(dat_turn$turn_e, dat_turn$turn_b, paired=TRUE)

#nestedness
nest_e <- designdist(eDNA_mat, "((pmax(b,c)-pmin(b,c)) / (a+b+c)) * (a / (a + (2 * pmin(b,c))))", abcd = T)
nest_b <- designdist(BS_mat, "((pmax(b,c)-pmin(b,c)) / (a+b+c)) * (a / (a + (2 * pmin(b,c))))", abcd = T)

n1 <- as.data.frame(matrix(nest_e, dimnames=list(t(outer(colnames(turn_e), rownames(turn_e), FUN=paste)), NULL)))
n2 <- as.data.frame(matrix(nest_b, dimnames=list(t(outer(colnames(turn_b), rownames(turn_b), FUN=paste)), NULL)))

dat_nest <- tibble(nest_e = n1, nest_b = n2)

plot_nest <- ggplot(data = dat_nest, aes(x = nest_b$V1, y = nest_e$V1)) +
  geom_point(size = 0.3) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1)) + 
  scale_y_continuous(limits = c(0, 1)) +
  annotate("segment", x = 0, xend = 0.9, y = 0, yend = 0.9) +
  annotate("text", x = 0.95, y = 0.95, label = "1:1") +
  ylab("Pairwise nestedness - eDNA") +
  xlab("Pairwise nestedness - beach seining") +
  annotate("text", x = 0.01, y = 1, label = "B")
plot_nest
ggsave("./Figures_tables/spatial_nest.png", 
       plot = plot_nest,
       width = 3, height = 3, units = "in")
#Jaccards
jac_e <- vegdist(eDNA_mat, "jaccard")
jac_b <- vegdist(BS_mat, "jaccard")

j1 <- as.data.frame(matrix(jac_e, dimnames=list(t(outer(colnames(turn_e), rownames(turn_e), FUN=paste)), NULL)))
j2 <- as.data.frame(matrix(jac_b, dimnames=list(t(outer(colnames(turn_b), rownames(turn_b), FUN=paste)), NULL)))

dat_jac <- tibble(jac_e = j1, jac_b = j2)

plot_jac <- ggplot(data = dat_jac, aes(x = jac_b$V1, y = jac_e$V1)) +
  geom_point(size = 0.3) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 1)) + 
  scale_y_continuous(limits = c(0, 1)) +
  annotate("segment", x = 0, xend = 0.9, y = 0, yend = 0.9) +
  annotate("text", x = 0.95, y = 0.95, label = "1:1") +
  ylab("Pairwise Jaccard's - eDNA") +
  xlab("Pairwise Jaccard's - beach seining") +
  annotate("text", x = 0.01, y = 1, label = "C")
plot_jac
ggsave("./Figures_tables/spatial_jac.png", 
       plot = plot_jac,
       width = 3, height = 3, units = "in")

#mantel test
mantel.rtest(turn_b, turn_e, nrepet = 999)
mantel.rtest(nest_b, nest_e, nrepet = 999)
mantel.rtest(jac_b, jac_e, nrepet = 999)


#PCoA ####

region <- data.frame(site = c("chp", "ssp", "wfb", "hab2", 
                              "hdo", "ppo", "pba", "hab4",
                              "fan1", "fan3", 
                              "gog1", "gog4",
                              "kis1", "kis2",
                              "sni1", "sni2",
                              "koe1", "koe3"),
                     region = c("1","1","1","1",
                                "2","2","2","2",
                                "3","3",
                                "3","3",
                                "3","3",
                                "4","4",
                                "4","4")) 
region2 <- data.frame(site = c("chp", "ssp", "wfb", "hab2", 
                              "hdo", "ppo", "pba", "hab4",
                              "fan1", "fan3", 
                              "gog1", "gog4",
                              "kis1", "kis2",
                              "sni1", "sni2",
                              "koe1", "koe3"),
                     region2 = c("1","1","1","1",
                                "2","2","2","2",
                                "1","1",
                                "1","1",
                                "1","1",
                                "2","2",
                                "2","2")) 

env1 <- env %>%
  merge(., region, by = "site") %>%
  merge(., region2, by = "site") %>%
  column_to_rownames(var = "dat_site") %>%
  mutate(year = factor(year, 
                       levels = c(2018,2019,2020), 
                       labels = c(2018,2019,2020))) %>%
  mutate(site_fac = factor(site, 
                         levels = c("chp", "ssp", "wfb", "hab2", 
                                    "hdo", "ppo", "pba", "hab4",
                                    "fan1", "fan3", "gog1", "gog4",
                                    "kis1", "kis2",
                                    "sni1", "sni2",
                                    "koe1", "koe3"), 
                         labels = c("chp", "ssp", "wfb", "hab2", 
                                    "hdo", "ppo", "pba", "hab4",
                                    "fan1", "fan3", "gog1", "gog4",
                                    "kis1", "kis2",
                                    "sni1", "sni2",
                                    "koe1", "koe3"))) %>%
  mutate(veg_05 = ifelse(subtidal_primary_cover > 5, "veg", "no_veg"))%>%
  mutate(veg_10 = ifelse(subtidal_primary_cover > 10, "veg", "no_veg"))%>%
  mutate(veg_15 = ifelse(subtidal_primary_cover > 15, "veg", "no_veg"))%>%
  mutate(veg_20 = ifelse(subtidal_primary_cover > 20, "veg", "no_veg"))
  
#PCoA beach seine
pcoa_b <- pcoa(turn_b)
pcoa_b
axes_b <- pcoa_b$vectors %>%
  as.data.frame() %>%
  merge(., env1, by = "row.names") %>%
  arrange(site_fac) %>%
  merge(.,aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~site_fac,.,mean),by="site_fac")
labels_b <- axes_b[c("site", "mean.x", "mean.y")] %>%
  distinct()
barplot(pcoa_b$values$Relative_eig[1:10])

d1 <- wascores(axes_b[c(3:4)], BS_mat) %>%
  as.data.frame() %>%
  mutate(vec_len = sqrt(Axis.1^2 + Axis.2^2))
biplot.pcoa(pcoa_b, BS_mat)


#BS eigenvalues - add tp plots below
pcoa_b$values$Relative_eig[1]
pcoa_b$values$Relative_eig[2]

#BS hulls
hull_b <- axes_b %>%
  group_by(region) %>%
  slice(chull(Axis.1, Axis.2))# %>%

#PCoA beach seine
pcoa_e <- pcoa(turn_e)
axes_e <- pcoa_e$vectors %>%
  as.data.frame() %>%
  merge(., env1, by = "row.names") %>%
  arrange(site_fac) %>%
  merge(.,aggregate(cbind(mean.x=Axis.1,mean.y=Axis.2)~site_fac,.,mean),by="site_fac")
labels_e <- axes_e[c("site", "mean.x", "mean.y")] %>%
  distinct()
barplot(pcoa_e$values$Relative_eig[1:10])

#eDNA eigenvalues
pcoa_e$values$Relative_eig[1]
pcoa_e$values$Relative_eig[2]

#eDNA hulls
hull_e <- axes_e %>%
  group_by(region) %>%
  slice(chull(Axis.1, Axis.2)) #%>%
#  filter(region == 1 | region == 2)

plot_e <- ggplot()+
  geom_polygon(data = hull_e,
               aes(x = Axis.1, y = Axis.2, fill = region),
               alpha = 0.3,
               show.legend = FALSE) +
  geom_point(data = axes_e, aes(x = Axis.1, y = Axis.2, shape = veg_15), size=1) +
  geom_point(data = axes_e, aes(x=mean.x,y=mean.y, shape = veg_15), size=0) +
  theme_classic()+
  geom_segment(data = axes_e, aes(x=mean.x, y=mean.y, xend=Axis.1, yend=Axis.2),linewidth = 0.1) +
  #geom_text(data = labels_e, aes(x=mean.x,y=mean.y, label = site), hjust = 0.5, vjust = -0.7) +
  # scale_shape_manual(values=c(15, 17, 18)) +
  xlab("PCoA-1 (40.1%)")+
  ylab("PCoA-2 (22.4%)")+ 
  theme(legend.position = "none",
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        plot.margin = margin(5,5,5,5))+ 
  #theme(legend.title=element_blank())+
  annotate("text", x = -0.264, y = 0.35, label = "eDNA", size = 3)

plot_e

plot_b <- ggplot() +
  geom_polygon(data = hull_b,
               aes(x = Axis.1, y = Axis.2, fill = region),
               alpha = 0.3,
               show.legend = FALSE) +
  geom_point(data = axes_b, aes(x = Axis.1, y = Axis.2, shape = veg_15), size=1) +
  geom_point(data = axes_b, aes(x=mean.x,y=mean.y, shape = veg_15), size=0) +
  theme_classic()+
  geom_segment(data = axes_b, aes(x=mean.x, y=mean.y, xend=Axis.1, yend=Axis.2),linewidth = 0.1) +
  #geom_text(data = labels_b, aes(x=mean.x,y=mean.y, label = site), hjust = 0.5, vjust = -0.7) +
  # scale_shape_manual(values=c(15, 17, 18)) +
  xlab("PCoA-1 (28.6%)")+
  ylab("PCoA-2 (24.7%)")+ 
  theme(legend.position = "bottom",
        legend.title=element_blank(),
        plot.margin = margin(0,5,5,5), 
        legend.margin = margin(-10,30,0,0),
        legend.text=element_text(size=8),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9),
        legend.text.align = 22,
       # legend.spacing.x = unit(0.1, "inches"),
        legend.box = "horizontal") + 
  scale_shape_discrete(breaks=c("no_veg", "veg"),
                       labels=c("vegetation cover < 15%", "vegetation cover > 15%"))+ 
  annotate("text", x = -0.29, y = 0.35, label = "beach seine", size = 3) 

plot_b

final <- 
ggarrange(plot_e, plot_b, 
          ncol = 1, nrow = 2,
          heights = c(0.92,1)) 

final

ggsave("./Figures_tables/final.png", 
       plot = final,
       width = 3.1, height = 4.3, units = "in")


#list of species inside and outside for eDNA and BS

r1 <- merge(long, env, by = "dat_site") %>%
  filter(site %in% c("chp", "wfb", "ssp", "hab2", "ppo", "pba", "hdo", "hab4")) %>%
  group_by(site, LCT_shared) %>%
  summarise(P_A_bs = sum(p_a_bs),
            P_A_eDNA = sum(p_a_eDNA)) %>%
  mutate(area = ifelse(site %in% c("chp", "wfb", "ssp", "hab2"), "outside", "inside"))
  #%>%
#  mutate(P_A_bs = if_else(P_A_bs == 0, 0, 1))
 
r2 <- r1 %>%
  select(-c(P_A_eDNA)) %>%
  rename(p_a = P_A_bs) %>%
  filter(p_a != 0)

in_bs <- r2 %>%
  filter(area == "inside") %>%
  ungroup() %>%
  distinct(LCT_shared)

out_bs <- r2 %>%
  filter(area == "outside") %>%
  ungroup() %>%
  distinct(LCT_shared)

r3 <- r1 %>%
  select(-c(P_A_bs)) %>%
  rename(p_a = P_A_eDNA) %>%
  filter(p_a != 0)

in_eDNA <- r3 %>%
  filter(area == "inside") %>%
  ungroup() %>%
  distinct(LCT_shared)

out_eDNA <- r3 %>%
  filter(area == "outside") %>%
  ungroup() %>%
  distinct(LCT_shared)

in_bs_unique <- in_bs$LCT_shared[!in_bs$LCT_shared %in% out_bs$LCT_shared]
out_bs_unique <- out_bs$LCT_shared[!out_bs$LCT_shared %in% in_bs$LCT_shared]

in_eDNA_unique <- in_eDNA$LCT_shared[!in_eDNA$LCT_shared %in% out_eDNA$LCT_shared]
out_eDNA_unique <- out_eDNA$LCT_shared[!out_eDNA$LCT_shared %in% in_eDNA$LCT_shared]

in_bs_unique
out_bs_unique
in_eDNA_unique
out_eDNA_unique

#GDMs ####

#filter for year. or not
env2 <- env1 %>%
  mutate(day_y_eDNA = lubridate::yday(date)) %>%
  mutate(day_y_BS = lubridate::yday(hakai_link_date)) %>%
  mutate(year_days = case_when(year == 2018 ~ 0,
                               year == 2019 ~ 365,
                               year == 2020 ~ 730)) %>%
  mutate(day_eDNA = day_y_eDNA + year_days) %>%
  mutate(day_BS = day_y_BS + year_days) %>%
  dplyr::select(-c("day_y_eDNA", "day_y_BS", "year_days"))

#turn this on to filter to a given year and month
#env2 <- env2 %>% 
# filter(year == 2018, month == 7)

eDNA_mat_gdm <- eDNA_mat %>% filter(row.names(eDNA_mat) %in% rownames(env2))
turn_e_gdm <- designdist(eDNA_mat_gdm, "2 * pmin(b,c) / (a + 2 * pmin(b,c))", abcd = T)

BS_mat_gdm <- BS_mat %>% filter(row.names(BS_mat) %in% rownames(env2))
turn_b_gdm <- designdist(BS_mat_gdm, "2 * pmin(b,c) / (a + 2 * pmin(b,c))", abcd = T)

#format distance matrix
mat_simple <- readRDS("Data/2022_10_31/derived_data/dist_mat2015.rds") #distance by water
mat_simple <- as.data.frame(as.matrix(mat_simple)) 
mat_simple <- mat_simple%>%
  filter(row.names(mat_simple) %in% env1$site) %>%
  dplyr::select(env1$site) %>%
  rownames_to_column(var = "site")

a1 <- env1[c(1)] %>%
  rownames_to_column("dat_site")
dist <- env1[c(1)] %>%
  rownames_to_column(var = "dat_site") %>%
  merge(., mat_simple, by = "site") %>%
  dplyr::select(-c("site")) %>%
  column_to_rownames(var = "dat_site") %>%
  t(.) %>%
  as.data.frame() %>%
  rownames_to_column(var = "site") %>%
  merge(a1,., by = "site") %>%
  dplyr::select(-c("site")) %>%
  rename("site" = "dat_site") %>%
  filter(site %in% rownames(env2)) %>%
  dplyr::select(c(1, rownames(env2)))
#dist[dist == 0] <- 5

t1 <- pivot_longer(dist, !site, names_to = "site_b", values_to = "dissim") %>%
  filter(dissim != 0) #%>%
#pivot_wider(names_from = site_b, values_from = dissim)

#summaries of distances
k1 <- dist %>%
  column_to_rownames("site") %>%
  as.matrix()
k2 <- as.data.frame(matrix(k1, dimnames=list(t(outer(colnames(k1), rownames(k1), FUN=paste)), NULL)))
hist(k2$V1, breaks = 50)

eDNA_GDM_mat <- 
  env2[c("lat", "long")] %>%
  merge(.,eDNA_mat, by = "row.names") %>%
  rename("site" = "Row.names")

BS_GDM_mat <- 
  env2[c("lat", "long")] %>%
  merge(.,BS_mat, by = "row.names") %>%
  rename("site" = "Row.names")

env3 <- env2 %>%
  rename("location" = "site") %>%
  rownames_to_column("site")

env_eDNA <- env3 %>%
  dplyr::select(c("site", "lat", "long", "exposure_num", "subtidal_primary_cover", "day_eDNA")) 
env_BS <- env3 %>%
  dplyr::select(c("site", "lat", "long", "exposure_num", "subtidal_primary_cover", "day_BS")) 

site_pair_e <- formatsitepair(eDNA_GDM_mat, bioFormat = 1, dist="bray", abundance=F, 
                            siteColumn= "site", XColumn = "long",YColumn = "lat", 
                            sppColumn=NULL, abundColumn=NULL, sppFilter=0, 
                            predData = env_eDNA, distPreds= list(as.matrix(dist)),
                            weightType="equal", custWeights=NULL, sampleSites=1) 
site_pair_e$distance <- as.vector(turn_e_gdm) #replace distance metric
#site_pair_e <- site_pair_e %>% #remove within site pairs
#  filter(s2.matrix_1 != 0)
site_pair_e <- site_pair_e %>%
  mutate(s2.matrix_1 = ifelse(s2.matrix_1 == 0, 5, s2.matrix_1))

# select number of splines for each predictor variable: only increasing distance i-splines improves deviance explained
splines <- c(3,3,3,4)
gdm_e <- gdm(site_pair_e, geo = F, splines = splines) # geo is false because I added a distance matrix (distance-by-water)
# view model result
str(gdm_e)
summary(gdm_e)
plot(gdm_e, plot.layout=c(3,2))

#this function (gdm.varImp) doesn't allow for differing spline numbers, here we use 4 (greatest from above)
#provides %deviance lost (if removed from model), p-values, model deviance explained, null deviance, GDM deviance, and intercept 
splines <- c(4,4,4,4)
set.seed(2)
var_imp_e <-  gdm.varImp(site_pair_e, geo = F, splines = splines, predSelect = F, 
                         nPerm = 500, pValue=0.05, parallel = T, cores = 7, sampleSites = 1, sampleSitePairs = 1,
                         outFile = NULL)
var_imp_e

#provides RMSE and observed-predicted correlation 
CV_e <- gdm.crossvalidation(spTable = site_pair_e, train.proportion=0.7, n.crossvalid.tests=20, geo=FALSE, splines=splines, knots=NULL)

#repeat for beach seining
site_pair_b <- formatsitepair(BS_GDM_mat, bioFormat = 1, dist="bray", abundance=F, 
                              siteColumn= "site", XColumn = "long",YColumn = "lat", 
                              sppColumn=NULL, abundColumn=NULL, sppFilter=0, 
                              predData = env_BS, distPreds= list(as.matrix(dist)),
                              weightType="equal", custWeights=NULL, sampleSites=1)
site_pair_b$distance <- as.vector(turn_b_gdm)
site_pair_b <- site_pair_b %>%
  mutate(s2.matrix_1 = ifelse(s2.matrix_1 == 0, 5, s2.matrix_1))

# select number of splines for each predictor variable, increasing any # of i-splines doesn't improve deviance explained by more than 1%
splines <- c(3,3,3,3)
gdm_b <- gdm(site_pair_b, geo = F, splines = splines) # geo is false because I added a distance matrix (distance-by-water)
# view model result
str(gdm_b)
summary(gdm_b)
plot(gdm_b, plot.layout=c(3,2))

var_imp_b <-  gdm.varImp(site_pair_b, geo = F, splines = splines, knots = NULL, predSelect = F, 
                         nPerm = 500, parallel = T, cores = 7, sampleSites = 1, sampleSitePairs = 1,
                         outFile = NULL)
var_imp_b
CV_b <- gdm.crossvalidation(spTable = site_pair_b, train.proportion=0.7, n.crossvalid.tests=20, geo=FALSE, splines=splines, knots=NULL)

#to get % deviance explained by each perdictor, run single variable models and take deviance explained from summary()

#get uncertainty estimates around i-splines
#do thes next 12 lines in order
uncertainly_e <- plotUncertainty(spTable = site_pair_e, sampleSites = 0.8, bsIters = 200, geo=FALSE, splines=c(3,3,3,4), knots=NULL,
                                 splineCol="black", errCol="grey80", plot.linewidth=2.0, plot.layout=c(3,2),
                                 parallel=T, cores=7, save = T)

error_e <- read.csv("gdm.plotUncertainy.csv") #written above


uncertainly_b <- plotUncertainty(spTable = site_pair_b, sampleSites = 0.8, bsIters = 200, geo=FALSE, splines=c(3,3,3,3), knots=NULL,
                                 splineCol="black", errCol="grey80", plot.linewidth=2.0, plot.layout=c(3,2),
                                 parallel=T, cores=7, save = T)

error_b <- read.csv("gdm.plotUncertainy.csv") #written above

#format error estimates for plotting
error_e <- error_e %>%
    select(c("exposure_num_fullModel_X", "exposure_num_minusSD_Y", "exposure_num_fullModel_Y", "exposure_num_plusSD_Y",
             "subtidal_primary_cover_fullModel_X", "subtidal_primary_cover_minusSD_Y", "subtidal_primary_cover_fullModel_Y", "subtidal_primary_cover_plusSD_Y",
             "day_eDNA_fullModel_X", "day_eDNA_minusSD_Y", "day_eDNA_fullModel_Y", "day_eDNA_plusSD_Y",
             "matrix_1_fullModel_X", "matrix_1_minusSD_Y", "matrix_1_fullModel_Y", "matrix_1_plusSD_Y")) %>%
  mutate(method = "eDNA")

error_e[error_e<0] <- 0
colnames(error_e) <- c("exp_x","exp_min","exp_y","exp_max",
                  "cover_x","cover_min","cover_y","cover_max",
                  "day_x","day_min","day_y","day_max",
                  "dist_x","dist_min","dist_y","dist_max", "method")

error_b <- error_b %>%
  select(c("exposure_num_fullModel_X", "exposure_num_minusSD_Y", "exposure_num_fullModel_Y", "exposure_num_plusSD_Y",
           "subtidal_primary_cover_fullModel_X", "subtidal_primary_cover_minusSD_Y", "subtidal_primary_cover_fullModel_Y", "subtidal_primary_cover_plusSD_Y",
           "day_BS_fullModel_X", "day_BS_minusSD_Y", "day_BS_fullModel_Y", "day_BS_plusSD_Y",
           "matrix_1_fullModel_X", "matrix_1_minusSD_Y", "matrix_1_fullModel_Y", "matrix_1_plusSD_Y")) %>%
  mutate(method = "beach_seine")
error_b[error_b<0] <- 0
colnames(error_b) <- c("exp_x","exp_min","exp_y","exp_max",
                       "cover_x","cover_min","cover_y","cover_max",
                       "day_x","day_min","day_y","day_max",
                       "dist_x","dist_min","dist_y","dist_max", "method")

error <- rbind(error_b, error_e)

exp_gdm <- 
  ggplot()  +
  geom_ribbon(data = error_e, aes(x = exp_x, ymin = exp_min, ymax = exp_max), fill = "#D55E00", alpha = 1) +
  geom_ribbon(data = error_b, aes(x = exp_x, ymin = exp_min, ymax = exp_max), fill = "#E69F00", alpha = 0.8) +
  geom_smooth(data = error_e, aes(x = exp_x, y = exp_y), method = lm, formula = y ~ splines::bs(x, 100), colour = "black", size = 0.5) +
  geom_smooth(data = error_b, aes(x = exp_x, y = exp_y), method = lm, formula = y ~ splines::bs(x, 100), colour = "black", size = 0.5) +
 theme_classic() +
  xlab("exposure (class)") +
  ylab("predicted turnover") +
  ylim(0,0.52) +
  theme(text = element_text(size = 12)) + 
  theme(legend.position = "none")
#  annotate("point", x = site_pair_b$s1.exposure_num, y = 0,colour = "black", size = 2.5, alpha=1, shape = "triangle")
exp_gdm
ggsave("./Figures_tables/GDM/exp_partial.png", 
       plot = exp_gdm,
       width = 5, height = 4, units = "in")
cover_gdm <- 
  ggplot()  +
  geom_ribbon(data = error_e, aes(x = cover_x, ymin = cover_min, ymax = cover_max), fill = "#D55E00", alpha = 1) +
  geom_ribbon(data = error_b, aes(x = cover_x, ymin = cover_min, ymax = cover_max), fill = "#E69F00", alpha = 0.8) +
  geom_smooth(data = error_e, aes(x = cover_x, y = cover_y), method = lm, formula = y ~ splines::bs(x, 100), colour = "black", size = 0.5) +
  geom_smooth(data = error_b, aes(x = cover_x, y = cover_y), method = lm, formula = y ~ splines::bs(x, 100), colour = "black", size = 0.5) +
  theme_classic() +
  xlab("vegetation cover (%)") +
  ylab("predicted turnover") +
  ylim(0,0.52) +
  theme(text = element_text(size = 12)) + 
  theme(legend.position = "none") 
cover_gdm
ggsave("./Figures_tables/GDM/cover_partial.png", 
       plot = cover_gdm,
       width = 5, height = 4, units = "in")

dist_gdm <- 
  ggplot()  +
  geom_ribbon(data = error_e, aes(x = dist_x, ymin = dist_min, ymax = dist_max), fill = "#D55E00", alpha = 1) +
  geom_ribbon(data = error_b, aes(x = dist_x, ymin = dist_min, ymax = dist_max), fill = "#E69F00", alpha = 0.8) +
  geom_smooth(data = error_e, aes(x = dist_x, y = dist_y), method = lm, formula = y ~ splines::bs(x, 100), colour = "black", size = 0.5) +
  geom_smooth(data = error_b, aes(x = dist_x, y = dist_y), method = lm, formula = y ~ splines::bs(x, 100), colour = "black", size = 0.5) +
  theme_classic() +
  xlab("distance (km)") +
  ylab("predicted turnover") +
  ylim(0,0.52) +
  theme(text = element_text(size = 12)) + 
  theme(legend.position = "none") +
  scale_x_continuous(expand = c(0, 0), breaks = c(10000, 20000,30000,40000,50000), labels = c("10", "20", "30", "40", "50"))#+
#  geom_point(data = site_pair_b, aes(x = s2.matrix_1, y = 0), shape = 3)
dist_gdm
ggsave("./Figures_tables/GDM/dist_partial.png", 
       plot = dist_gdm,
       width = 5, height = 4, units = "in")

day_gdm <- 
  ggplot()  +
  geom_ribbon(data = error_e, aes(x = day_x, ymin = day_min, ymax = day_max), fill = "#D55E00", alpha = 1) +
  geom_ribbon(data = error_b, aes(x = day_x, ymin = day_min, ymax = day_max), fill = "#E69F00", alpha = 0.8) +
  geom_smooth(data = error_e, aes(x = day_x, y = day_y), method = lm, formula = y ~ splines::bs(x, 100), colour = "black", size = 0.5) +
  geom_smooth(data = error_b, aes(x = day_x, y = day_y), method = lm, formula = y ~ splines::bs(x, 100), colour = "black", size = 0.5) +
  theme_classic() +
  xlab("date") +
  ylab("predicted turnover") +
  theme(text = element_text(size = 12)) + 
    theme(legend.position = "none")+
  ylim(0,0.52) +
  scale_x_continuous(breaks=c(200, 525, 850), labels = c("2018","2019", "2020")) +
  annotate("rect", xmin = 420, xmax = 550, ymin = 0.44, ymax = 0.5, alpha = 0.8,fill = "#E69F00")+
  annotate("rect", xmin = 420, xmax = 550, ymin = 0.36, ymax = 0.42, alpha = 0.8,fill = "#D55E00")+
  annotate("segment", x = 420, xend = 550, y = 0.47, yend = 0.47, linewidth = 0.5)+
  annotate("segment", x = 420, xend = 550, y = 0.39, yend = 0.39, linewidth = 0.5) +
  annotate("text", x = 745, y = 0.473, label = "beach seine model") +
  annotate("text", x = 692, y = 0.393, label = "eDNA model") +
  annotate("rect", xmin = 400, xmax = 940, ymin = 0.34, ymax = 0.52, colour = "gray", alpha = 0)

day_gdm
ggsave("./Figures_tables/GDM/day_partial.png", 
       plot = day_gdm,
       width = 5, height = 4, units = "in")

gdm_all <- 
ggarrange(dist_gdm, day_gdm, cover_gdm, exp_gdm,  
          ncol = 2, nrow = 2,
          heights = c(1,1)) 
gdm_all
ggsave("./Figures_tables/GDM/gdm_all.png", 
       plot = gdm_all,
       width = 8, height = 6, units = "in")

