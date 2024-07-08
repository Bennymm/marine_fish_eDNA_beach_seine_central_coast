#comparing alpha diversity between eDNA and BS 
#Model taxonomic and functional richness as functions of habitat heterogeneity and method
#Test correlation between richness from each method - see species compared (almost 
# always more with eDNA (~10 high variance))
# packages and data ####
library(tidyverse)
library(lubridate)
library(here)
library(lme4)
library(scales)
library(vegan)
library(AICcmodavg)
library(performance)
library(MASS)
library(MuMIn)

#long data to calculate richness
long <- read_rds("no_occupancy_model/Data/2022_10_31/shared_long.rds")
#envonmental data
env <- read_rds("no_occupancy_model/Data/2022_10_31/environmental.rds")

#calculate richness ####
rich <- long %>%
  dplyr::select(c("dat_site", "p_a_bs", "p_a_eDNA"))  %>%
  group_by(dat_site) %>%
  mutate(across(c("p_a_bs", "p_a_eDNA"), sum)) %>%
  `colnames<-` (c("dat_site","rich_b","rich_e")) %>%
  distinct() %>%
  mutate(rich_diff = rich_e - rich_b)


#plots of richness differences ####
richxrich <-
ggplot(rich, aes(x = rich_b, y = rich_e)) +
  geom_point() + 
  theme_classic() +
  xlab("LCT richness - beach seining") +
  ylab("LCT richness - eDNA") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 36)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 36)) +
  annotate("segment", x = 0, xend = 29, y = 0, yend = 29) +
  annotate("text", x = 30, y = 30, label = "1:1")
richxrich


rich_long <- rich %>%
  dplyr::select(-c("rich_diff")) %>%
  pivot_longer(!dat_site, names_to = "method", values_to = "richness")

richdiff <- 
ggplot(rich_long, aes(x = method, y = richness, group = dat_site, colour =)) +
  geom_point() +
  geom_line() +
  ylab("Taxonomic\nrichness") +
  xlab("Method") + 
  theme_classic()+ 
  scale_x_discrete(labels=c("rich_b" = "beach seine", "rich_e" = "eDNA"))
richdiff



richdiffdensity <- 
ggplot(rich, aes(x = rich_diff)) +
  geom_histogram(binwidth = 3, aes( y = ..density..), fill = "white", color = "black", size = 0.5) +
  geom_density() +
  theme_classic() +
  geom_vline(xintercept = 0) +
  xlab(" LCT richness difference (eDNA - beach seine)") +
  ylab("density of\nrichness differences") 
richdiffdensity

#summaries
mean(rich$rich_b)
sd(rich$rich_b)
mean(rich$rich_e)
sd(rich$rich_e)
mean(rich$rich_diff)
sd(rich$rich_diff)

# models of richness differences ####

data <- merge(rich, env, by = "dat_site")

library(nlme)
# correlation structures: https://norcalbiostat.github.io/AppliedStatistics_notes/specifying-correlation-structures.html
#                         https://www.rdocumentation.org/packages/nlme/versions/3.1-160/topics/corClasses


  
data1 <- data %>%
  mutate(site_num = as.numeric(factor(site))) 


#habitat proximity
mm01 <- gls(rich_diff ~  seagrass + kelp,
            correlation = corSymm(form = ~1|site),
            data = data1)
mm02 <- gls(rich_diff ~ rockyshore + water25m + slope_subtidal, 
            correlation = corSymm(form = ~1|site),
            data = data1)
mm03 <- gls(rich_diff ~ freshwater, 
            correlation = corSymm(form = ~1|site),
            data = data1)
#habitat diversity
mm04 <- gls(rich_diff ~ h1000m + h100m, 
            correlation = corSymm(form = ~1|site),
            data = data1)
#seawater turnover
mm05 <- gls(rich_diff ~ silt_percent,
            correlation = corSymm(form = ~1|site),
            data = data1)
#temporal offset
mm06 <- gls(rich_diff ~ dat_diff,
            correlation = corSymm(form = ~1|site),
            data = data1)
#mixed
mm07 <- gls(rich_diff ~ h100m + h1000m +
              silt_percent, 
            correlation = corSymm(form = ~1|site),
            data = data1)

vif(mm01)
vif(mm02)
vif(mm03)
vif(mm04)
vif(mm05)
vif(mm06)
vif(mm07)

options(scipen=999)
model_sel <- as.data.frame(model.sel(mm01, mm02, mm03, mm04, mm05, mm06, mm07)) %>%
  rownames_to_column(var = "model")
model_sel
write_csv(model_sel, "no_occupancy_model/Figures_tables/richmodels_summary.csv")
options(scipen=0)
model_performance(mm07)
summary(mm07)

vif(mm07)

newdata <- data1[c("silt_percent", "h1000m", "h100m")]
predict <- predict(mm07, newdata = newdata)
plot(x = data1$rich_diff, y = predict)



#standardized beta coefficients
data2 <- scale(data1[c("h100m", "h1000m", "silt_percent")]) %>%
  as.data.frame() %>%
  cbind(data1[c("rich_diff", "site")],.)
mm07_std <- gls(rich_diff ~ h100m + h1000m +
                  silt_percent, 
                correlation = corSymm(form = ~1|site),
                data = data2)
coefficients(mm07_std)

confint <- as.data.frame(confint(mm07_std))
ceof <- as.data.frame(mm07_std$coefficients)
coef_confint <- cbind(ceof, confint)
coef_confint <- rownames_to_column(coef_confint, "vars")
colnames(coef_confint) <- c("var", "beta", "lo", "hi")

ggplot(data=coef_confint,aes(x=var, y=beta, ymin = lo, ymax = hi)) +
  geom_point() +
  geom_errorbar(width = 0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), 
        axis.title.x=element_blank()) +
  ylab("Standardized\ncoefficient") +
  geom_hline(yintercept = 0) +
  scale_x_discrete(labels = c("intercept", "# of habitats\nwithin 1000m", "# of habitats\nwithin 100m", "% fine\n sediment"))


#correlation between habitat richness and habitat richness
data3 <- data1 %>%
  group_by(site, h100m, h1000m, silt_percent)%>%
  summarise(rich_bs = mean(rich_b),
            rich_e = mean(rich_e))

cor.test(x = data3$h100m, y= data3$rich_bs, method = "spearman")
cor.test(x = data3$h100m, y= data3$rich_e, method = "spearman")
cor.test(x = data3$silt_percent, y= data3$rich_bs, method = "spearman")
cor.test(x = data3$silt_percent, y= data3$rich_e, method = "spearman")
cor.test(x = data3$h1000m, y= data3$rich_bs, method = "spearman")
cor.test(x = data3$h1000m, y= data3$rich_e, method = "spearman")



#plot of survey design
m1 <- env %>%
  dplyr::select(c("date", "dat_site", "hakai_link_date")) %>%
  pivot_longer(!dat_site, names_to = "date_type", values_to = "date") %>%
  merge(., env[c("site", "dat_site")], by = "dat_site")

ggplot(data = m1, aes(x = date, y = site, group = dat_site)) +
  geom_point() +
  geom_line() +
  theme_classic()
