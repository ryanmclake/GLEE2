library(dplyr)
library(vroom)
library(readr)
library(ggplot2)
library(patchwork)
library(sf)
library(units)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(tidyverse)
library(sp)
library(rworldmap)
library(arrow)
library(nimble)
library(coda)

d <- vroom::vroom("./data/organized_data_to_append/global_lake_res_CH4_emission_DB_all.csv")

length(unique(d$primary_data_source))

duplicates <- d[duplicated(d$primary_data_source) | duplicated(d$primary_data_source, fromLast = TRUE), ]

length(unique(duplicates$primary_data_source))

repeats <- unique(duplicates$primary_data_source)

data_wo_repeats <- d %>% filter(!primary_data_source %in% repeats)
data_w_repeats <- d %>% filter(primary_data_source %in% repeats)

johnson_res <- data_w_repeats %>% filter(published_data_product == "Johnson et al 2021") %>% select(lat, lon, ch4_diff, ch4_ebu, published_data_product, primary_data_source, obs_years)


johnson_lake <- data_w_repeats %>% filter(published_data_product == "Johnson et al 2022") %>% select(lat, lon, ch4_diff, ch4_ebu, published_data_product, primary_data_source, obs_years)


prairie_res <- data_w_repeats %>% filter(published_data_product == "Prairie et al 2021") %>% select(lat, lon, ch4_diff, ch4_ebu, published_data_product, primary_data_source, obs_years)


rosentreter_both <- data_w_repeats %>% filter(published_data_product == "Rosentreter et al 2020") %>% select(lat, lon, ch4_diff, ch4_ebu, published_data_product, primary_data_source, obs_years) 


kuhn_lake <- data_w_repeats %>% filter(published_data_product == "Kuhn et al 2021") %>% select(lat, lon, ch4_diff, ch4_ebu, published_data_product, primary_data_source, obs_years)



# Link Johnson to Prairie

link1 <- johnson_res %>%
  cbind(prairie_res[st_nearest_feature(johnson_res, prairie_res),]) %>%
  mutate(dist = (st_distance(geometry, geometry.1, by_element = T))*0.001) %>%
  st_drop_geometry(.) %>% 
  arrange(dist) %>%
  filter(ref == ref.1)

compare_Jr_and_Pr_diff <- link1 %>% select(ch4_diff, ch4_diff.1, published_data_product.1, primary_data_source.1) %>% 
  mutate(ch4_diff.1 = ifelse(ch4_diff.1 == 0, NA, ch4_diff.1)) %>%
  na.omit(.) %>%
  ggplot(., aes(x = ch4_diff.1, y = ch4_diff, color = primary_data_source.1))+
  geom_point(alpha = 1, size = 4, pch = 19)+
  theme_classic()+
  #geom_smooth(method = lm, se = T, alpha = 0.1)+
  geom_abline(slope = 1)+
  ylab("Johnson et al 2021 - Diffusion")+
  xlab("Prairie et al 2018 - Diffusion")+
  labs(title = "Johnson vs. Prairie - Diffusion")

compare_Jr_and_Pr_ebu <- link1 %>% select(ch4_ebu, ch4_ebu.1, published_data_product.1, primary_data_source.1) %>% 
  mutate(ch4_ebu.1 = ifelse(ch4_ebu.1 == 0, NA, ch4_ebu.1)) %>%
  na.omit(.) %>%
  ggplot(., aes(x = ch4_ebu.1, y = ch4_ebu, color = primary_data_source.1))+
  geom_point(alpha = 1, size = 4, pch = 19)+
  theme_classic()+
  #geom_smooth(method = lm, se = T, alpha = 0.1)+
  geom_abline(slope = 1)+
  ylab("Johnson et al 2021 - Ebullition")+
  xlab("Prairie et al 2018 - Ebullition")+
  labs(title = "Johnson vs. Prairie - Ebullition")

# Link Johnson_Res to Rosentreter All

link2 <- johnson_res %>% left_join(., rosentreter_both, by = c("lat","obs_years")) %>% na.omit()

compare_Jr_and_R_diff <- link2 %>% select(ch4_diff.x, ch4_diff.y, published_data_product.x, primary_data_source.x) %>% 
  na.omit(.) %>%
  ggplot(., aes(x = ch4_diff.x, y = ch4_diff.y, color = primary_data_source.x))+
  geom_point(alpha = 1, size = 4, pch = 19)+
  theme_classic()+
  #geom_smooth(method = lm, se = T, alpha = 0.1)+
  geom_abline(slope = 1)+
  ylab("Johnson et al 2021 - Diffusion")+
  xlab("Rosentreter et al 2020 - Diffusion")+
  labs(title = "Johnson vs. Rosentreter - Diffusion")

compare_J_and_R_ebu <- link2 %>% select(ch4_ebu, data_source, ch4_ebu.1, data_source.1, ref) %>% 
  mutate(ch4_diff.1 = ifelse(ch4_ebu.1 == 0, NA, ch4_ebu.1)) %>%
  na.omit(.) %>%
  ggplot(., aes(x = ch4_ebu.1, y = ch4_ebu, color = ref))+
  geom_point(alpha = 1, size = 4, pch = 19)+
  theme_classic()+
  #geom_smooth(method = lm, se = T, alpha = 0.1)+
  geom_abline(slope = 1)+
  ylab("Johnson et al 2021 - Ebullition")+
  xlab("Rosentreter et al 2020 - Ebullition")+
  labs(title = "Johnson vs. Rosentreter - Ebullition")

# Link Johnson to Kuhn

link3 <- johnson_lake %>% left_join(., kuhn_lake, by = c("lat","obs_years")) %>% na.omit()


compare_J_and_K_diff <- link3 %>% select(ch4_diff.x, ch4_diff.y, published_data_product.x, primary_data_source.x) %>% 
  na.omit(.) %>%
  ggplot(., aes(x = ch4_diff.x, y = ch4_diff.y, color = primary_data_source.x))+
  geom_point(alpha = 1, size = 4, pch = 19)+
  theme_classic()+
  #geom_smooth(method = lm, se = T, alpha = 0.1)+
  geom_abline(slope = 1)+
  ylab("Johnson et al 2021 - Diffusion")+
  xlab("Kuhn et al 2020 - Diffusion")+
  labs(title = "Johnson vs. Kuhn - Diffusion")

compare_J_and_K_ebu <- link3 %>% select(ch4_ebu, data_source, ch4_ebu.1, data_source.1, ref) %>% 
  mutate(ch4_diff.1 = ifelse(ch4_ebu.1 == 0, NA, ch4_ebu.1)) %>%
  na.omit(.) %>%
  ggplot(., aes(x = ch4_ebu.1, y = ch4_ebu, color = ref))+
  geom_point(alpha = 1, size = 4, pch = 19)+
  theme_classic()+
  #geom_smooth(method = lm, se = T, alpha = 0.1)+
  geom_abline(slope = 1)+
  ylab("Johnson et al 2021 - Ebullition")+
  xlab("Kuhn et al 2020 - Ebullition")+
  labs(title = "Johnson vs. Kuhn - Ebullition")

# Link Rosentreter to Kuhn

link4 <- rosentreter %>%
  cbind(kuhn[st_nearest_feature(rosentreter, kuhn),]) %>%
  mutate(dist = (st_distance(geometry, geometry.1, by_element = T))*0.001) %>%
  st_drop_geometry(.) %>%
  arrange(dist) %>%
  filter(ref == ref.1)

compare_R_and_K_diff <- link4 %>% select(ch4_diff, data_source, ch4_diff.1, data_source.1, ref) %>% 
  mutate(ch4_diff.1 = ifelse(ch4_diff.1 == 0, NA, ch4_diff.1)) %>%
  ggplot(., aes(x = ch4_diff.1, y = ch4_diff, color = ref))+
  geom_point(alpha = 1, size = 4, pch = 19)+
  theme_classic()+
  #geom_smooth(method = lm, se = T, alpha = 0.1)+
  geom_abline(slope = 1)+
  ylab("Rosentreter et al 2021 - Diffusion")+
  xlab("Kuhn et al 2020 - Diffusion")+
  labs(title = "Rosentreter vs. Kuhn - Diffusion")

compare_R_and_K_ebu <- link4 %>% select(ch4_ebu, data_source, ch4_ebu.1, data_source.1, ref) %>% 
  mutate(ch4_diff.1 = ifelse(ch4_ebu.1 == 0, NA, ch4_ebu.1)) %>%
  ggplot(., aes(x = ch4_ebu.1, y = ch4_ebu, color = ref))+
  geom_point(alpha = 1, size = 4, pch = 19)+
  theme_classic()+
  #geom_smooth(method = lm, se = T, alpha = 0.1)+
  geom_abline(slope = 1)+
  ylab("Rosentreter et al 2021 - Ebullition")+
  xlab("Kuhn et al 2020 - Ebullition")+
  labs(title = "Rosentreter vs. Kuhn - Ebullition")

# Link Rosentreter to Prairie

link5 <- rosentreter %>%
  cbind(prairie[st_nearest_feature(rosentreter, prairie),]) %>%
  mutate(dist = (st_distance(geometry, geometry.1, by_element = T))*0.001) %>%
  st_drop_geometry(.) %>%
  arrange(dist) %>%
  filter(ref == ref.1)

compare_R_and_P_diff <- link5 %>% select(ch4_diff, data_source, ch4_diff.1, data_source.1, ref) %>% 
  mutate(ch4_diff.1 = ifelse(ch4_diff.1 == 0, NA, ch4_diff.1)) %>%
  ggplot(., aes(x = ch4_diff.1, y = ch4_diff, color = ref))+
  geom_point(alpha = 1, size = 4, pch = 19)+
  theme_classic()+
  #geom_smooth(method = lm, se = T, alpha = 0.1)+
  geom_abline(slope = 1)+
  ylab("Rosentreter et al 2021 - Diffusion")+
  xlab("Prairie et al 2020 - Diffusion")+
  labs(title = "Rosentreter vs.Prairie - Diffusion")








