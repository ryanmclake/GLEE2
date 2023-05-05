library(dplyr)
library(vroom)
library(readr)
library(ggplot2)
library(patchwork)
library(sf)
library(units)

# Load in the dataset
d1 <- vroom::vroom("./data/organized_data_to_append/global_lake_res_DB_countries.csv")

# determine how many unique refs there are
length(unique(d1$ref))

# group the refs within data sources and get their global mean
d3 <- d1 %>% group_by(ref, lat, lon, obs_year, country, continent, num_months_sampled) %>%
  summarize(mean_ebu = mean(ch4_ebu, na.rm = T),
            mean_diff = mean(ch4_diff, na.rm = T),
            sd_ebu = sd(ch4_ebu, na.rm = T),
            sd_diff = sd(ch4_diff, na.rm = T),
            mean_obs_wtemp_k = mean(effective_obs_wtemp_k, na.rm = T),
            sd_obs_wtemp_k = sd(effective_obs_wtemp_k, na.rm = T),
            mean_cum_radiance = mean(cumulative_radiance, na.rm = T),
            sd_cum_radiance = sd(cumulative_radiance, na.rm = T),
            mean_TP_ugL = mean(TP_ugL, na.rm = T),
            sd_TP_ugL = sd(TP_ugL, na.rm = T),
            mean_littoral_area_km2 = mean(littoral_area, na.rm = T),
            sd_littoral_area_km2 = sd(littoral_area, na.rm = T),
            mean_surf_area_km2 = mean(surf_area_k, na.rm = T),
            sd_surf_area_km2 = sd(surf_area_k, na.rm = T),
            mean_depth_m = mean(mean_depth, na.rm = T),
            sd_mean_depth_m = sd(mean_depth, na.rm = T),
            max_depth_m = mean(max_depth, na.rm = T),
            sd_max_depth_m = sd(max_depth, na.rm = T))



month_calc <- d3 %>% group_by(num_months_sampled) %>%
  summarize(sd_ebu = sd(mean_ebu, na.rm = T),
            sd_diff = sd(mean_diff, na.rm = T))


# group the number of months and get their global mean and SD
d4 <- d3 %>% group_by(num_months) %>%
  summarize(sd_ebu = sd(mean_ebu, na.rm = T),
            sd_diff = sd(mean_diff, na.rm = T))

# confirm that ref length is the same as original data file
length(unique(d3$ref))

duplicates <- d3[duplicated(d3$ref) | duplicated(d3$ref, fromLast = TRUE), ]

duplicates_J_and_R <- duplicates %>% filter(data_source %in% c("Johnson et al 2021","Rosentreter et al 2020"))
J_R_repeats <- c(unique(duplicates_J_and_R$ref))

# plot all of the repeated data points and their corresponding value based on 
# the original data source
duplicates_ebu <- duplicates %>% 
  ggplot(., aes(x=ref, y=mean_ebu, group = data_source, fill = data_source))+
  geom_point(size = 4, pch = 21, color = "black", alpha = 0.5)+
  theme_classic()+
  ylab("Ebullition Rate (mgCH4/m2*d)")+
  xlab("")+
  labs(title = paste0("N repeat = ",length(unique(duplicates$ref))))+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

# plot all of the repeated data points and their corresponding value based on 
# the original data source
duplicates_diff <- duplicates %>%
  ggplot(., aes(x=ref, y=mean_diff, group = data_source, fill = data_source))+
  geom_point(size = 4, pch = 21, color = "black", alpha = 0.5)+
  theme_classic()+
  ylab("Diffusion Rate (mgCH4/m2*d)")+
  xlab("")+
  labs(title = paste0("N repeat = ",length(unique(duplicates$ref))))+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

# combine the ebullition and diffusion plots
repeat_compare <- duplicates_diff/duplicates_ebu
repeat_compare

# get the mean of all the duplicated vales across the three data sources
repeats <- unique(duplicates$ref)

# filter all of the repeats from the original data source and take summary stats
# DO THIS INSTEAD OF TAKING MEANS FROM d2. This way we don't violate Jensen's 
# inequality

# Make two data products, one with repeats and one without
data_wo_repeats <- d1 %>% filter(!ref %in% repeats)
data_w_repeats <- d1 %>% filter(ref %in% repeats)

johnson <- data_w_repeats %>% filter(data_source == "Johnson et al 2021") %>% select(lat, lon, ch4_diff, ch4_ebu, ref, data_source) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")

prairie <- data_w_repeats %>% filter(data_source == "Prairie et al 2018") %>% select(lat, lon, ch4_diff, ch4_ebu, ref, data_source) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")

rosentreter <- data_w_repeats %>% filter(data_source == "Rosentreter et al 2020") %>% select(lat, lon, ch4_diff, ch4_ebu, ref, data_source) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")




# Link Johnson to Prairie - only one DP (Tremblay et al 2005)
link1 <- johnson %>%
  cbind(prairie[st_nearest_feature(johnson, prairie),]) %>%
  mutate(dist = (st_distance(geometry, geometry.1, by_element = T))*0.001) %>%
  st_drop_geometry(.) %>%
  arrange(dist) %>%
  filter(ref == ref.1) %>%
  mutate(mean_diff_repeat = (ch4_diff+ch4_diff.1)/2)

compare_J_and_P_diff <- ggplot(link1, aes(x = ch4_diff.1, y = ch4_diff, fill = ref))+
  geom_point(color = "black", alpha = 1, size = 4, pch = 21)+
  theme_classic()+
  #geom_smooth(method = lm, se = T, alpha = 0.1)+
  geom_abline(slope = 1)+
  ylab("Johnson et al 2021 - Diffusion")+
  xlab("Prairie et al 2018 - Diffusion")+
  labs(title = "Johnson vs. Prairie - Diffusion")


# Link Johnson to Rosentreter - only one DP (Tremblay et al 2005)
link2 <- johnson %>%
  cbind(rosentreter[st_nearest_feature(johnson, rosentreter),]) %>%
  mutate(dist = (st_distance(geometry, geometry.1, by_element = T))*0.001) %>%
  st_drop_geometry(.) %>%
  arrange(dist) %>%
  filter(ref == ref.1) %>%
  mutate(mean_diff_repeat = (ch4_diff+ch4_diff.1)/2,
         mean_ebu_repeat = (ch4_ebu+ch4_ebu.1)/2)

compare_J_and_R_diff <- link2 %>% select(ch4_diff, data_source, ch4_diff.1, data_source.1, ref) %>% 
  mutate(ch4_diff.1 = ifelse(ch4_diff.1 == 0, NA, ch4_diff.1)) %>%
  na.omit(.) %>%
  ggplot(., aes(x = ch4_diff.1, y = ch4_diff, color = ref))+
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




# Link Prairie to Rosentreter - only one DP (Tremblay et al 2005)
link3 <- prairie %>%
  cbind(rosentreter[st_nearest_feature(prairie, rosentreter),]) %>%
  mutate(dist = (st_distance(geometry, geometry.1, by_element = T))*0.001) %>%
  st_drop_geometry(.) %>%
  arrange(dist) %>%
  filter(ref == ref.1)%>%
  mutate(mean_diff_repeat = (ch4_diff+ch4_diff.1)/2,
         mean_ebu_repeat = (ch4_ebu+ch4_ebu.1)/2)

compare_P_and_R_diff <- link3 %>% select(ch4_diff, data_source, ch4_diff.1, data_source.1, ref) %>% 
  mutate(ch4_diff.1 = ifelse(ch4_diff.1 == 0, NA, ch4_diff.1)) %>%
  na.omit(.) %>%
  ggplot(., aes(x = ch4_diff.1, y = ch4_diff, color = ref))+
  geom_point(alpha = 1, size = 4, pch = 19)+
  theme_classic()+
  #geom_smooth(method = lm, se = T, alpha = 0.1)+
  geom_abline(slope = 1)+
  ylab("Prairie et al 2021 - Diffusion")+
  xlab("Rosentreter et al 2020 - Diffusion")+
  labs(title = "Prairie vs. Rosentreter - Diffusion")

compare_P_and_R_ebu <- link3 %>% select(ch4_ebu, data_source, ch4_ebu.1, data_source.1, ref) %>% 
  mutate(ch4_diff.1 = ifelse(ch4_ebu.1 == 0, NA, ch4_ebu.1)) %>%
  na.omit(.) %>%
  ggplot(., aes(x = ch4_ebu.1, y = ch4_ebu, color = ref))+
  geom_point(alpha = 1, size = 4, pch = 19)+
  theme_classic()+
  #geom_smooth(method = lm, se = T, alpha = 0.1)+
  geom_abline(slope = 1)+
  ylab("Prairie et al 2021 - Ebullition")+
  xlab("Rosentreter et al 2020 - Ebullition")+
  labs(title = "Prairie vs. Rosentreter - Ebullition")

linked_compare <- (compare_J_and_P_diff + plot_spacer())/(compare_J_and_R_diff+compare_J_and_R_ebu)/(compare_P_and_R_diff+compare_P_and_R_ebu)
linked_compare


