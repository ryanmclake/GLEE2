
rm(list=ls())
gc()

### NOTE ### JAGS needs to be downloaded in order to run this analysis beginning to end
# Download the JAGS application at this link:
# https://sourceforge.net/projects/mcmc-jags/files/
# CMD (or Control) + Click to follow link


# install packages needed for the analysis

if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, MCMCvis, lubridate, tidybayes, readr,
               ncdf4, reshape2, zoo, patchwork, hydroGOF, viridis,
               imputeTS, devtools, scales, forecast, coda, rjags, 
               R2jags, gridExtra, maps, hexbin, rnaturalearth, sf, nlstools)


# read in data set that has been linked to GLCP Hydrobasin climate data
base <- read_csv("./source_data/GLEE_data_with_GLCP_HWSD_link.csv")


filtered_lakes <- base %>%
  ## --> quantify the number of total samples from each site by multiplying the
  ## --> number of months sampled by the number of sites collecting CH4 emissions
  mutate(tot_sampling_events = num_months_sampled * num_sites_sampled) %>%
  group_by(geometry, waterbody_type) %>%
  ## --> calculate a waterbody_id column for each individual lakes
  ## --> this will convey observation uncertainty by individual lakes
  mutate(waterbody_id = cur_group_id()) %>%
  ## --> remove any lakes that have only one measurement
  #filter(any(length(waterbody_id)>1)) %>%
  ungroup() %>%
  ## --> assign the temperature used for the model to find either the observed air temp 
  ## --> or the air temperature as determined by the GCLP hydrobasin
  mutate(temp_for_model_K = ifelse(is.na(mean_temp_k), effective_obs_wtemp_k, mean_temp_k)) %>%
  group_by(waterbody_id) %>%
  mutate(temp_for_model_K = ifelse(n_distinct(effective_obs_wtemp_k) == 0, mean_temp_k, effective_obs_wtemp_k))

# Write the output of this lake thinning to a new file
write_csv(filtered_lakes, "./output/filtered_GLEE_w_GLCP_link.csv")

# Calculate the standard deviation of the ebullition observations of each lake 
error_ebu <- filtered_lakes %>% select(waterbody_id, waterbody_type,num_months_sampled, 
                                       num_sites_sampled, tot_sampling_events, ch4_ebu,
                                       temp_for_model_K, surf_area_k, T_OC, mean_sw_wm2) %>% 
    ## --> filter only the ebullition observations
    ## --> Group by waterbody ID
    group_by(waterbody_id) %>%
    ## --> calculate the SD based on the number of total months sampled, the number
    ## --> of sites sampled, and the product of sites and months sampled
    mutate(sd_time = (sqrt(sum(abs(ch4_ebu - mean(ch4_ebu))^2))/num_months_sampled),
           sd_space = (sqrt(sum(abs(ch4_ebu - mean(ch4_ebu))^2))/num_sites_sampled),
           temp_for_model_C = temp_for_model_K - 273.15,
           sd_time = ifelse(sd_time == 0, ch4_ebu*3, sd_time),
           sd_space = ifelse(sd_space == 0, ch4_ebu*3, sd_space))

# Write it to a file so you don't always have to rerun this script if R Fails
write_csv(error_ebu, "./output/filtered_GLEE_ebullition_w_GLCP_link_and_ERROR.csv")

# Calculate the standard deviation of the ebullition observations of each lake 
error_diff <- filtered_lakes %>% select(waterbody_id, waterbody_type, num_months_sampled, 
                                        num_sites_sampled, tot_sampling_events, ch4_diff,
                                        temp_for_model_K, surf_area_k, T_OC, mean_sw_wm2) %>% 
  ## --> filter only the diffusion observations
  ## --> Group by waterbody ID
  group_by(waterbody_id) %>%
  ## --> calculate the SD based on the number of total months sampled, the number
  ## --> of sites sampled, and the product of sites and months sampled
  mutate(sd_time = (sqrt(sum(abs(ch4_diff - mean(ch4_diff))^2))/num_months_sampled),
         sd_space = (sqrt(sum(abs(ch4_diff - mean(ch4_diff))^2))/num_sites_sampled),
         temp_for_model_C = temp_for_model_K - 273.15,
         sd_time = ifelse(sd_time == 0, ch4_diff*2, sd_time),
         sd_space = ifelse(sd_space == 0, ch4_diff*2, sd_space))

# Write it to a file so you don't always have to rerun this script if R Fails
write_csv(error_diff, "./output/filtered_GLEE_diffusion_w_GLCP_link_and_ERROR.csv")

# Make a global plot of the site locations used in the analysis

base2 <- read_csv("./source_data/GLEE_data_with_GLCP_link.csv")

## --> Download the globe from R Natural Earth
world <- ne_countries(scale = "medium", returnclass = "sf")

## --> Assign the flux pathway to the lakes and reservoirs
k <- base2 %>%
  select(lat, lon, waterbody_type, ch4_ebu, ch4_diff) %>%
  group_by(lat, lon, waterbody_type) %>%
  summarize_all(funs(mean)) %>%
  mutate(`Emission` = ifelse(!is.na(ch4_ebu),"Both Fluxes", "Only Diffusion"),
         `Emission` = ifelse(is.na(ch4_diff),"Only Ebullition", `Emission`)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(waterbody_type == "lake" | waterbody_type == "reservoir") %>%
  rename(`Waterbody Type` = waterbody_type)

## --> Make a global plot
ggplot() +
  geom_sf(data = world, lwd = 0.3, color = "black", fill = "white")+
  xlab("Longitude") + ylab("Latitude") +
  labs(title = paste0("  ",length(k$`Waterbody Type`)," Lakes and Reservoirs"))+
  geom_sf(data = k, size = 2.5, aes(shape = Emission, color = `Waterbody Type`))+
  scale_color_manual(values = c("purple", "cyan"))+
  theme_void()+
  theme(legend.position = c(0.11, 0.5),
        legend.direction = "vertical",
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = element_text(size=9),
        legend.key.height  = unit(.5, 'cm'),
        legend.key.width =  unit(.3, 'cm'))

k 
