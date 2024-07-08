#test statistics for Gamma-level errors in Beach seine and eDNA surveys

library(tidyverse)
library(nlme)
library(car)
library(scales)
library(performance)
data_raw <- read_csv("./Figures_tables/LCT_table_habitat_preferences.csv")
data1 <- filter(data_raw, gamma_error != "eDNA only") %>%
  mutate(eDNA_pa = ifelse(prev_eDNA > 0, 1, 0))

# eDNA detection and biomass density
ggplot(data = data1, aes(x = max_biodensBS, y = eDNA_pa)) +
  geom_point() +
  scale_x_continuous(trans = pseudo_log_trans(base = 2), 
                     breaks = c(0,1,10,100,1000,10000,100000),
                     labels = c("0","1   ", "10 ", expression(10^{"2"}), expression(10^{"3"}), expression(10^{"4"}), expression(10^{"5"})),
                     limits = c(0, 1000))

bs_model = glm(eDNA_pa ~ max_biodensBS, family = "binomial", data = data1)
summary(bs_model)
performance(bs_model)
AIC(bs_model)
anova(bs_model, test = "Chisq")
Anova(bs_model, type = 3)

# eDNA only and habitat
#data
data2 <- data_raw %>%
  filter(expected != "na") %>%
  select(c("gamma_error", "expected")) %>%
  mutate(gamma_error = ifelse(gamma_error == "eDNA only", "eDNA only", "other"))

table(data2$gamma_error, data2$expected)

#binomial
a1 <- data2 %>%
  group_by(gamma_error, expected) %>%
  summarise(length = length(expected))
a1
48/(48+7)
2/26
binom.test(2,26,0.87)

