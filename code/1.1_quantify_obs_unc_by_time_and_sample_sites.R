if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, MCMCvis, lubridate, tidybayes,
               ncdf4, reshape2, zoo, patchwork, hydroGOF, viridis,
               imputeTS, devtools, scales, forecast, coda, rjags, R2jags,gridExtra,
               maps, hexbin, rnaturalearth, sf)

base <- read_csv("./data/organized_data_to_append/GLEE_data_with_GLCP_link.csv") %>%
  mutate(tot_sampling_events = num_months_sampled * num_sites_sampled) %>%
  group_by(lat, lon) %>%
  mutate(waterbody_id = cur_group_id()) %>%
  filter(any(length(waterbody_id)>1)) %>%
  ungroup() %>%
  mutate(temp_for_model_K = ifelse(is.na(mean_temp_k), effective_obs_wtemp_k, mean_temp_k),
         temp_for_model_C = temp_for_model_K-273.15)

error_ebu <- base %>% select(waterbody_id, waterbody_type, num_months_sampled, num_sites_sampled, tot_sampling_events, ch4_ebu, temp_for_model_C) %>%
    group_by(waterbody_id) %>%
    mutate(sd_time = (sqrt(sum(abs(ch4_ebu - mean(ch4_ebu))^2))/num_months_sampled),
           sd_space = (sqrt(sum(abs(ch4_ebu - mean(ch4_ebu))^2))/num_sites_sampled),
           sd_space_time = (sqrt(sum(abs(ch4_ebu - mean(ch4_ebu))^2))/tot_sampling_events)) %>%
  filter(sd_space_time > 0) %>%
  na.omit(.)

error_diff <- base %>% select(waterbody_id, waterbody_type, num_months_sampled, num_sites_sampled, tot_sampling_events, ch4_diff, temp_for_model_C) %>%
  group_by(waterbody_id) %>%
  mutate(sd_time = (sqrt(sum(abs(ch4_diff - mean(ch4_diff))^2))/num_months_sampled),
         sd_space = (sqrt(sum(abs(ch4_diff - mean(ch4_diff))^2))/num_sites_sampled),
         sd_space_time = (sqrt(sum(abs(ch4_diff - mean(ch4_diff))^2))/tot_sampling_events)) %>%
  filter(sd_space_time > 0) %>%
  na.omit(.)









world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

k <- base %>%
  select(lat, lon, ch4_ebu, ch4_diff) %>%
  group_by(lat, lon) %>%
  summarize_all(funs(mean)) %>%
  mutate(emission = ifelse(!is.na(ch4_ebu),"Both", "Diffusion"),
         emission = ifelse(is.na(ch4_diff),"Ebullition", emission)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") 

ggplot() +
  geom_sf(data = world, lwd = 0.3, color = "black", fill = "white")+
  xlab("Longitude") + ylab("Latitude") +
  labs(title = paste0("    521 Sites"))+
  geom_sf(data = k, size = 2.5, aes(shape = emission, color = emission))+
  scale_color_viridis_d(option = "D")+
  theme_void()+
  theme(legend.position = c(0.11, 0.5),
        legend.direction = "vertical",
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = element_text(size=9),
        legend.key.height  = unit(.5, 'cm'),
        legend.key.width =  unit(.3, 'cm'))

