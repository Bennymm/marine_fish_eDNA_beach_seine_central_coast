library(tidyverse)

long <- read_rds("Data/2022_10_31/shared_long.rds")
#envonmental data
env <- read_rds("Data/2022_10_31/environmental.rds")


long1 <- long %>%
  merge(., env[c("dat_site","h100m", "site", "year", "month")], by = "dat_site") #%>%
  filter(year == 2018)


t1 <- long1 %>%
  group_by(dat_site) %>%
  mutate(rank_bio = order(order(grams_per_m2, decreasing=TRUE)),
         rank_abund = order(order(abundance, decreasing=TRUE))) %>%
  filter(grams_per_m2 > 0.000000000000001) %>%
  mutate(
         p_a_eDNA = as.character(p_a_eDNA))

t2 <- t1 %>%
  filter(p_a_eDNA == 0)

ggplot(t1, aes(x = rank_bio, y = grams_per_m2, group = dat_site, colour = h100m, shape = p_a_eDNA)) +
  #  geom_point(size = 2) +
  geom_line()+
  geom_point(data = t2, aes(x = rank_bio, y = grams_per_m2, fill = "black"))+
  scale_y_log10() +
  theme_classic()


pd <- position_dodge(0)
ggplot() +
  #  geom_point(size = 2) +
  geom_line(data = t1, aes(x = rank_bio, y = grams_per_m2, group = dat_site, colour = h100m),position = pd)+
  scale_colour_gradient(low = "red", high = "green")+
  geom_point(data = t2, aes(x = rank_bio, y = grams_per_m2, shape = p_a_eDNA),position = pd)+
  scale_shape_manual(values=c(1, 3))+
  scale_y_continuous(trans = 'log10',limits = c(0.0001,1e3), 
                                         breaks = c(0.001,0.01,0.1,1,10,100,1000),
                                          labels = c(expression(10^{"-3"}),
                                                     expression(10^{"-2"}),
                                                     expression(10^{"-1"}),
                                                     "1", 
                                                     "10",
                                                     expression(10^{"2"}), 
                                                     expression(10^{"3"}))) +
  theme_classic()+ 
  theme(legend.position = c(0.8,0.7))  + 
  ylab(bquote('Biomass density '(g / m^2))) +
  xlab("Rank abundance")






h1 <- filter(t1, h100m == 1)
h2 <- filter(t1, h100m == 2)
h3 <- filter(t1, h100m == 3)
h4 <- filter(t1, h100m == 4)

ggplot() +
  geom_point(data = h1, aes(x = rank_bio, y = grams_per_m2, colour = "blue"))+
  geom_smooth(data = h1, aes(x = rank_bio, y = grams_per_m2, colour = "blue"),method = lm, formula = y ~ splines::bs(x, 1), se = FALSE) +

  geom_point(data = h2, aes(x = rank_bio, y = grams_per_m2, colour = "green"))+
  geom_smooth(data = h2, aes(x = rank_bio, y = grams_per_m2, colour = "green"),method = lm, formula = y ~ splines::bs(x, 1), se = FALSE) +
 
  geom_point(data = h3, aes(x = rank_bio, y = grams_per_m2, colour = "black"))+
  geom_smooth(data = h3, aes(x = rank_bio, y = grams_per_m2, colour = "black"),method = lm, formula = y ~ splines::bs(x, 1), se = FALSE) +
  
  geom_point(data = h4, aes(x = rank_bio, y = grams_per_m2, colour = "yellow"))+
  geom_smooth(data = h4, aes(x = rank_bio, y = grams_per_m2, colour = "yellow"),method = lm, formula = y ~ splines::bs(x, 1), se = FALSE) +
  
  scale_y_log10() +
  theme_classic()
    
ggplot()+
  
  scale_y_log10() +
  theme_classic() +
  geom_point(data = h1, aes(x = rank_bio, y = grams_per_m2, colour = "blue", size = 0.5))+
  geom_smooth(data = h1, aes(x = rank_bio, y = grams_per_m2, colour = "blue"),method = lm, formula = y ~ splines::bs(x, 3), se = FALSE) +
  
  geom_point(data = h4, aes(x = rank_bio, y = grams_per_m2, colour = "green", size = 0.5)) +
  geom_smooth(data = h4, aes(x = rank_bio, y = grams_per_m2, colour = "green"),method = lm, formula = y ~ splines::bs(x, 3), se = FALSE) 




#black point around species not seen with eDNA
#do this with subset of species that eDNA did see

#plot richness in each method against richness at both 1000m and 100m

#calculate richness ####
rich <- long %>%
  dplyr::select(c("dat_site", "p_a_bs", "p_a_eDNA"))  %>%
  group_by(dat_site) %>%
  mutate(across(c("p_a_bs", "p_a_eDNA"), sum)) %>%
  `colnames<-` (c("dat_site","rich_b","rich_e")) %>%
  distinct() %>%
  mutate(rich_diff = rich_e - rich_b) %>%
  merge(., env[c("dat_site","h100m", "h1000m")], by = "dat_site")

rich1 <- rich %>% 
  pivot_longer(cols = c("h100m", "h1000m"), names_to = "scale", values_to = "hab_rich")

ggplot(rich1, aes(x = hab_rich, y = rich_b, group = scale, colour = scale, size = 1)) +
  geom_jitter(width = 0.05, height = 2) +
  theme_classic()

ggplot(rich1, aes(x = hab_rich, y = rich_e, group = scale, colour = scale, size = 1)) +
  geom_jitter(width = 0.05, height = 2)+
  theme_classic()


rich2 <- rich %>%
  pivot_longer(cols = c("rich_e", "rich_b"), names_to = "method", values_to = "richness")%>%
  pivot_longer(cols = c("h100m", "h1000m"), names_to = "scale", values_to = "hab_rich")

h100m <- rich2 %>% filter(scale == "h100m")
h1000m <- rich2 %>% filter(scale == "h1000m")


h100m_bs <- h100m %>% filter(method == "rich_b")
h100m_bs <- h100m %>% filter(method == "rich_b")

#lm_bs <- lm(richness ~ hab_rich, data = h100m_bs)
#summary(lm_bs)
#cor.test(x = h100m_bs$hab_rich, y= h100m_bs$richness, method = "spearman")
#plot(x = h100m_bs$hab_rich, y= h100m_bs$richness)


pd <- position_dodge(0.5)
ggplot(h100m, aes(x = hab_rich, y = richness, group = dat_site, shape = method, colour = dat_site)) +
  geom_point(position = pd, size = 3)+
  geom_line(position = pd) +
  theme_classic()+
  geom_smooth(aes(group = method, colour = "black"),
              method = "nls", formula = y ~ a * x + b, se = F,
              method.args = list(start = list(a = 0.1, b = 0.7)))


ggplot(h1000m, aes(x = hab_rich, y = richness, group = dat_site, shape = method, colour = dat_site)) +
  geom_point(position = pd, size = 3)+
  geom_line(position = pd) +
  theme_classic()+
  geom_smooth(aes(group = method, colour = "black"),
              method = "nls", formula = y ~ a * x + b, se = F,
              method.args = list(start = list(a = 0.1, b = 0.7)))


















