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

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

k <- read_csv("./data/organized_data_to_append/global_lake_res_DB_combined_refs_WWF_add.csv") %>%
  select(lat, lon)

f <- read_csv("./data/organized_data_to_append/global_lake_res_DB_combined_refs_WWF_add.csv") %>%
  select(-geometry) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")

month_error <- f %>% group_by(num_months_sampled) %>%
  summarise(ebu_sd = sd(mean_ebu, na.rm = T),
            diff_sd = sd(mean_diff, na.rm = T),
            temp_sd = sd(mean_obs_wtemp_k, na.rm = T))

h <- f %>%
  mutate(ebu_sd = ifelse(mean_ebu > 0, 1, sd_ebu)) %>%
  mutate(ebu_sd = ifelse(mean_ebu > 2, 1.5, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(mean_ebu > 5, 3, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(mean_ebu > 10, 5, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(mean_ebu > 50, 10, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(mean_ebu > 100, 40, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(mean_ebu > 200, 80, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(mean_ebu > 500, 200, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(mean_ebu > 1000, 500, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(mean_ebu > 3000, 1600, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(mean_ebu > 5000, 2500, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(mean_ebu > 7000, 3000, ebu_sd)) %>%
  mutate(diff_sd = ifelse(mean_diff > 0, 0.5, sd_diff)) %>%
  mutate(diff_sd = ifelse(mean_diff > 2, 0.8, diff_sd)) %>%
  mutate(diff_sd = ifelse(mean_diff > 5, 2, diff_sd)) %>%
  mutate(diff_sd = ifelse(mean_diff > 10, 3, diff_sd)) %>%
  mutate(diff_sd = ifelse(mean_diff > 50, 10, diff_sd)) %>%
  mutate(diff_sd = ifelse(mean_diff > 100, 30, diff_sd)) %>%
  mutate(diff_sd = ifelse(mean_diff > 200, 60, diff_sd)) %>%
  mutate(diff_sd = ifelse(mean_diff > 500, 100, diff_sd)) %>%
  mutate(diff_sd = ifelse(mean_diff > 1000, 300, diff_sd)) %>%
  mutate(diff_sd = ifelse(mean_diff > 3000, 500, diff_sd)) %>%
  mutate(diff_sd = ifelse(mean_diff > 5000, 700, diff_sd)) %>%
  select(-sd_ebu, -sd_diff)
  
bind_cols(h, k) %>% write_csv(., "./data/organized_data_to_append/global_lake_res_DB_combined_refs_WWF_sd_added.csv")


obs_diffusive_flux <- h %>% select(mean_diff, diff_sd, waterbody_type) %>% na.omit(.) %>%
  ggplot() +
  geom_sf(data = world, lwd = 0.3, color = "black", fill = "white")+
  xlab("Longitude") + ylab("Latitude") +
  labs(title = paste0("    N diffusive observations = 1166"))+
  geom_sf(data = h, pch = 21, color = "black", aes(size = mean_diff, fill = diff_sd))+
  scale_size("Methane Diffusion Rate")+
  scale_fill_gradient2(midpoint=0, low="tan1", mid="dodgerblue4",
                       high="red", space ="Lab", na.value="grey",
                       name = "Standard Deviation") +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = c(0.11, 0.5),
        legend.direction = "vertical",
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = element_text(size=9),
        legend.key.height  = unit(.5, 'cm'),
        legend.key.width =  unit(.3, 'cm'))

obs_ebullition_flux <- h %>% select(mean_ebu, ebu_sd) %>% na.omit(.) %>%
  ggplot() +
  geom_sf(data = world, lwd = 0.3, color = "black", fill = "white")+
  xlab("Longitude") + ylab("Latitude") +
  labs(title = paste0("    N ebullition observations = 611"))+
  geom_sf(data = h, pch = 21, color = "black", aes(size = mean_ebu, fill = ebu_sd))+
  scale_size("Methane Ebullition Rate")+
  scale_fill_gradient2(midpoint=0, low="tan1", mid="dodgerblue4",
                       high="red", space ="Lab", na.value="grey",
                       name = "Standard Deviation") +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = c(0.11, 0.5),
        legend.direction = "vertical",
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = element_text(size=9),
        legend.key.height  = unit(.5, 'cm'),
        legend.key.width =  unit(.3, 'cm'))
