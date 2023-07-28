
source("./code/ebullition_functions.R")

library(dplyr)
library(readr)

d <- vroom::vroom(paste0("/Users/ryanmcclure/Documents/GLEE2.1/data/countries/Japan.csv"), col_names = F, delim = ",") %>%
  rename(year = X1, month= X2,
         hylak_id = X3, lat = X4,
         lon = X5, continent = X6, country = X7,
         total_km2 = X8, waterbody_type = X9,
         mean_temp_k = X10, ice_cover_mean = X11) %>%
  mutate(ice_cover_mean = ifelse(is.na(ice_cover_mean),0,ice_cover_mean)) %>%
  na.omit(.) %>%
  select(-continent) %>%
  group_by(hylak_id, lat, lon, total_km2) %>%
  mutate(FOA_mean_ebullition = ifelse(waterbody_type == 1, mapply(FOA_mean_lake, x=mean_temp_k-273.15), mapply(FOA_mean_res, x=mean_temp_k-273.15)),
         FOA_time_low_ebu = ifelse(waterbody_type == 1, mapply(FOA_time_low_lake, x=mean_temp_k-273.15), mapply(FOA_time_low_res, x=mean_temp_k-273.15)),
         FOA_time_high_ebu = ifelse(waterbody_type == 1, mapply(FOA_time_high_lake, x=mean_temp_k-273.15), mapply(FOA_time_high_res, x=mean_temp_k-273.15)),
         FOA_space_low_ebu = ifelse(waterbody_type == 1, mapply(FOA_space_low_lake, x=mean_temp_k-273.15), mapply(FOA_space_low_res, x=mean_temp_k-273.15)),
         FOA_space_high_ebu = ifelse(waterbody_type == 1, mapply(FOA_space_high_lake, x=mean_temp_k-273.15), mapply(FOA_space_high_res, x=mean_temp_k-273.15)),
         FOA_mean_error_low_ebu = ifelse(waterbody_type == 1, mapply(FOA_mean_error_low_lake, x=mean_temp_k-273.15), mapply(FOA_mean_error_low_res, x=mean_temp_k-273.15)),
         FOA_mean_error_high_ebu = ifelse(waterbody_type == 1, mapply(FOA_mean_error_high_lake, x=mean_temp_k-273.15), mapply(FOA_mean_error_high_res, x=mean_temp_k-273.15)),
         FOA_coefficient_low_ebu = ifelse(waterbody_type == 1, mapply(FOA_coefficient_low_lake, x=mean_temp_k-273.15), mapply(FOA_coefficient_low_res, x=mean_temp_k-273.15)),
         FOA_coefficient_high_ebu = ifelse(waterbody_type == 1, mapply(FOA_coefficient_high_lake, x=mean_temp_k-273.15), mapply(FOA_coefficient_high_res, x=mean_temp_k-273.15)),
         FOA_mean_ebullition = ifelse(ice_cover_mean > 50, 0, FOA_mean_ebullition),
         FOA_time_low_ebu = ifelse(ice_cover_mean > 50, 0, FOA_time_low_ebu),
         FOA_time_high_ebu = ifelse(ice_cover_mean > 50, 0, FOA_time_high_ebu),
         FOA_space_low_ebu = ifelse(ice_cover_mean > 50, 0, FOA_space_low_ebu),
         FOA_space_high_ebu = ifelse(ice_cover_mean > 50, 0, FOA_space_high_ebu),
         FOA_mean_error_low_ebu = ifelse(ice_cover_mean > 50, 0, FOA_mean_error_low_ebu),
         FOA_mean_error_high_ebu = ifelse(ice_cover_mean > 50, 0, FOA_mean_error_high_ebu),
         FOA_coefficient_low_ebu = ifelse(ice_cover_mean > 50, 0, FOA_coefficient_low_ebu),
         FOA_coefficient_high_ebu = ifelse(ice_cover_mean > 50, 0, FOA_coefficient_high_ebu),
         FOA_mean_ebullition = ifelse(FOA_mean_ebullition < 0, 0, FOA_mean_ebullition),
         FOA_time_low_ebu = ifelse(FOA_time_low_ebu < 0, 0, FOA_time_low_ebu),
         FOA_time_high_ebu = ifelse(FOA_time_high_ebu < 0, 0, FOA_time_high_ebu),
         FOA_space_low_ebu = ifelse(FOA_space_low_ebu < 0, 0, FOA_space_low_ebu),
         FOA_space_high_ebu = ifelse(FOA_space_high_ebu < 0, 0, FOA_space_high_ebu),
         FOA_mean_error_low_ebu = ifelse(FOA_mean_error_low_ebu < 0, 0, FOA_mean_error_low_ebu),
         FOA_mean_error_high_ebu = ifelse(FOA_mean_error_high_ebu < 0, 0, FOA_mean_error_high_ebu),
         FOA_coefficient_low_ebu = ifelse(FOA_coefficient_low_ebu < 0, 0, FOA_coefficient_low_ebu),
         FOA_coefficient_high_ebu = ifelse(FOA_coefficient_high_ebu < 0, 0, FOA_coefficient_high_ebu),
         'date' = lubridate::make_date(year = year, month = month)) %>%
  ungroup(.) %>%
  select(-year, -month) %>%
  readr::write_csv(., "/Users/ryanmcclure/Documents/GLEE2.1/output/J_ebullition_variability.csv")



