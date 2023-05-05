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


d2 %>% write_csv("./data/organized_data_to_append/global_lake_res_DB_countries.csv")

# Load in the dataset
d_diff <- d3 %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  select(mean_diff, sd_diff, geometry)

d_ebu <- d3 %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  select(mean_ebu, sd_ebu, geometry) %>%
  mutate(sd_ebu = ifelse(is.na(sd_ebu), 10, sd_ebu))

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

obs_diffusive_flux <- 
  ggplot() +
  geom_sf(data = world, lwd = 0.3, color = "black", fill = "white")+
  xlab("Longitude") + ylab("Latitude") +
  labs(title = paste0("    N diffusive observations = ", length(d_diff$mean_diff)))+
  geom_sf(data = d_diff, pch = 21, color = "black", aes(size = mean_diff, fill = sd_diff))+
  scale_fill_gradient2(midpoint=0, low="tan1", mid="dodgerblue4",
                       high="red", space ="Lab", na.value="grey",
                       name = "**Diffusive Methane Emissions** <br>√(mgCH4/m2*d)") +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = c(0.11, 0.5),
        legend.direction = "vertical",
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = element_text(size=9),
        legend.key.height  = unit(.5, 'cm'),
        legend.key.width =  unit(.3, 'cm'))

obs_ebullition_flux <- 
  ggplot() +
  geom_sf(data = world, lwd = 0.3, color = "black", fill = "white")+
  xlab("Longitude") + ylab("Latitude") +
  labs(title = paste0("    N ebullition observations = ", length(d_ebu$mean_ebu)))+
  geom_sf(data = d_ebu, pch = 21, color = "black", aes(size = mean_ebu, fill = sd_ebu))+
  scale_fill_gradient2(midpoint=0, low="tan1", mid="dodgerblue4",
                       high="red", space ="Lab", na.value="grey",
                       name = "**Ebullitive Methane Emissions** <br>√(mgCH4/m2*d)") +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = c(0.11, 0.5),
        legend.direction = "vertical",
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = element_text(size=9),
        legend.key.height  = unit(.5, 'cm'),
        legend.key.width =  unit(.3, 'cm'))
