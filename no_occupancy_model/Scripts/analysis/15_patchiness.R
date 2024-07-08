#vegetation and exposure patchiness

library(tidyverse)
library(here)
library(stringr)

#vegetation patch sizes - select largest dimension
seagrass_dims <- read_csv("./Data/2022_10_31/spatial/seagrass_dims.csv")
macro_dims <- read_csv("./Data/2022_10_31/spatial/macro_dims.csv")
sea1 <- seagrass_dims %>%
  mutate(max_dim = if_else(width > height, width, height))
mac1 <- macro_dims %>%
  mutate(max_dim = if_else(width > height, width, height))
hist(sea1$max_dim, breaks = 40)
hist(mac1$max_dim, breaks = 40)


#exposure patch size 
  #
exp <- read_csv("./Data/2022_10_31/spatial/trial.csv")
#format geometry
e1 <- exp %>%
  dplyr::select(c("WKT", "EXP_FINAL", "FEAT_LEN")) %>%
  distinct()%>%
  mutate(geometry = str_remove_all(WKT, "[MULTILINESTRING()]")) %>%
  mutate(dd = as.list(str_split(geometry,",")))
e2 <- e1$dd

#find first and last lat/long for each segment of coastline
first <- c()
last <- c()
for (i in 1:length(e2)) {
  a <- first(e2[[i]])
  b <- last(e2[[i]])
  first <- c(first, a)
  last <- c(last, b)
  }
e1$first = first
e1$last = last

#format first and last
fun1 <- function(x) (as.numeric())
fun2 <- function(x) (round(x, digits = 0))
e3 <- e1 %>%
  mutate(first = str_remove(first, "[ ]")) %>%
#  mutate(last = str_remove(last, "[ ]")) %>%
  separate_wider_delim(first, " ", names = c("f_x", "f_y")) %>%
  separate_wider_delim(last, " ", names = c("l_x", "l_y")) %>%
  mutate(f_x = as.numeric(f_x),
         f_y = as.numeric(f_y),
         l_x = as.numeric(l_x),
         l_y = as.numeric(l_y),
         f_x = round(f_x),
         f_y = round(f_y),
         l_x = round(l_x),
         l_y = round(l_y),
         first = paste(f_x, f_y, sep = "_"),
         last = paste(l_x, l_y, sep = "_")) %>%
  select(c("EXP_FINAL", "FEAT_LEN", "first", "last")) %>%
  filter(!is.na(EXP_FINAL))


firstandlast <- c(e3$first, e3$last) %>%
  as.data.frame() %>%
  distinct()

#make a for loop that goes through each first and last in e3, merging together if they have the same exposure
exp_vec <- c("VP", "P", "SP", "SE", "E", "VE")
df <- data.frame()
for (j in exp_vec) {
  P1 <- e3 %>% filter(EXP_FINAL == j)
  P2 <- c(P1$first, P1$last) %>%
    as.data.frame() %>%
    distinct() 
    for (i in 1:nrow(P2)){
  
      m1 <- filter(P1, first == P2[i,1] | last == P2[i,1])
      if (nrow(m1) == 2) {
        m2 <- data.frame(EXP_FINAL = j, FEAT_LEN = sum(m1$FEAT_LEN), first = m1$first[2], last = m1$last[1])
        P3 <- filter(P1, first != P2[i,1] & last != P2[i,1])
        P1 <- rbind(m2,P3)
      }
    }
  df <- rbind(df,P1)
}


#format and merge
exp <- df %>% mutate(feature = "exposure") %>%
  dplyr::select("feature", "FEAT_LEN") %>%
  `colnames<-`(c("feature", "length"))
zos <- s2 %>% mutate(feature = "seagrass") %>%
  dplyr::select("feature", "max_dim") %>%
  `colnames<-`(c("feature", "length"))
macro <- m2 %>% mutate(feature = "macrocystis") %>%
  dplyr::select("feature", "max_dim") %>%
  `colnames<-`(c("feature", "length"))
all <- rbind(exp, zos, macro) %>%
  filter(length > 10) #filter out few very small patches

#summaries
all_summary <- all %>% group_by(feature) %>% summarise(mean = median(length),
                                                       perc_10th = quantile(length, 0.1),
                                                       perc_90th = quantile(length, 0.9))

#plotting
ggplot(all, aes(x = length, colour = feature, fill = feature))  +
  theme_classic() +
  geom_histogram(aes(y=..density..),position = "identity", alpha = 0.0, bins = 50) +
  geom_density(alpha = 0.2) +
  scale_x_continuous(trans = "log", breaks = c(1,10,100,1000,10000)) +
  xlab("length (log10 scale)") +
  scale_colour_discrete(labels=c('exposure (n=9609)', 'macrocystis (n=3089)', 'seagrass (n=6881)')) +
  scale_fill_discrete(labels=c('exposure (n=9609)', 'macrocystis (n=3089)', 'seagrass (n=6881)')) + 
  theme(legend.position = c(0.8, 0.8))

ggplot(all, aes(x = length, colour = feature, fill = feature))  +
  theme_classic() +
  geom_density(alpha = 0.2) +
  scale_x_continuous(trans = "log", breaks = c(1,10,100,1000,10000)) +
  xlab("length (meters; log10 scale)") +
  scale_colour_discrete(labels=c('exposure (n=9609)', 'macrocystis (n=3089)', 'seagrass (n=6881)')) +
  scale_fill_discrete(labels=c('exposure (n=9609)', 'macrocystis (n=3089)', 'seagrass (n=6881)')) + 
  theme(legend.position = c(0.8, 0.8))









