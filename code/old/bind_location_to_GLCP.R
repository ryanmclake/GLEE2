

dat <- read_csv("./data/organized_data_to_append/global_lake_res_DB_with_WWF.csv")

country <- dat %>% filter(country == "Sweden") %>% select(-geometry) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  rename(year = obs_year)

glcp <- read_csv_arrow(paste0("./data/glcp_countries/sweden.csv"),   
               quote = "\"",
               escape_double = TRUE,
               escape_backslash = FALSE,
               schema = NULL,
               col_names = F,
               col_types = NULL,
               col_select = NULL,
               na = c("", "NA"),
               quoted_na = TRUE,
               skip_empty_rows = TRUE,
               skip = 0L,
               parse_options = NULL,
               convert_options = NULL,
               read_options = NULL,
               as_data_frame = TRUE,
               timestamp_parsers = NULL) %>% rename(year = f0,
                                                    month=f1,
                                                    hylak_id=f2,
                                                    centr_lat=f3,
                                                    centr_lon=f4,
                                                    continent=f5,
                                                    country=f6,
                                                    bsn_lvl=f7,
                                                    hybas_id=f8,
                                                    mean_monthly_precip_mm=f9,
                                                    total_precip_mm=f10,
                                                    mean_annual_temp_k=f11,
                                                    pop_sum=f12,
                                                    seasonal_km2=f13,
                                                    permanent_km2=f14,
                                                    total_km2=f15,
                                                    lake_name=f16,
                                                    lake_type=f17,
                                                    lake_area=f18,
                                                    shore_dev=f19,
                                                    vol_total=f20,
                                                    vol_res=f21,
                                                    vol_src=f22,
                                                    depth_avg=f23,
                                                    res_time=f24,
                                                    elevation=f25,
                                                    slope_100=f26,
                                                    wshd_area=f27,
                                                    pour_long=f28,
                                                    pour_lat=f29,
                                                    sub_area=f30,
                                                    mean_spec_humidity=f31,
                                                    mean_precip_mm=f32,
                                                    sum_precip_mm=f33,
                                                    mean_temp_k=f34,
                                                    mean_totcloud_pct=f35,
                                                    mean_sw_wm2=f36,
                                                    mean_lw_wm2=f37,
                                                    above_ratio_cutoff=f38,
                                                    ice_cover_min=f39,
                                                    ice_cover_max=f40,
                                                    ice_cover_mean=f41,
                                                    ice_cover_median=f42,
                                                    ice_cover_binary_min=f43,
                                                    ice_cover_binary_max=f44,
                                                    ice_cover_binary_mean=f45,
                                                    ice_cover_binary_median=f46,
                                                    ice_cover_count=f47,
                                                    snow_km2=f48) %>%
  group_by(year, hylak_id, centr_lat, centr_lon,continent, country, bsn_lvl, hybas_id, 
           elevation, sub_area, lake_area, slope_100, depth_avg, shore_dev) %>%
  #summarizing observations at each lake by year 
  #using 'median' as our summary statistic 
  summarize(mean_monthly_precip_mm = mean(mean_monthly_precip_mm, na.rm = T),
            total_precip_mm = mean(total_precip_mm, na.rm = T),
            mean_annual_temp_k = mean(mean_annual_temp_k, na.rm = T),
            pop_sum = mean(pop_sum, na.rm = T),
            seasonal_km2 = mean(seasonal_km2, na.rm = T),
            permanent_km2 = mean(permanent_km2, na.rm = T),
            total_km2 = mean(total_km2, na.rm = T),
            ice_cover_min = mean(ice_cover_min, na.rm = T),
            ice_cover_max = mean(ice_cover_max, na.rm = T),
            ice_cover_mean = mean(ice_cover_mean, na.rm = T),
            ice_cover_mean = mean(ice_cover_mean, na.rm = T),
            ice_cover_count = mean(ice_cover_count, na.rm = T),
            mean_spec_humidity = mean(mean_spec_humidity, na.rm = T),
            mean_sw_wm2 = mean(mean_sw_wm2, na.rm = T),
            mean_lw_wm2 = mean(mean_lw_wm2, na.rm = T)) %>%
  #renaming columns back to sensible names
  #collecting so it is in arrow table format
  collect() %>% st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")

test <- do.call('rbind', lapply(split(country, 1:nrow(country)), function(x) {
  st_join(x, glcp[glcp$year == unique(x$year),], join = st_nearest_feature)
}))


