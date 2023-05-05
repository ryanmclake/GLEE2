
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "data")
Dir.Shapes <- file.path(Dir.Data, "shapes")

if (!file.exists(file.path(Dir.Shapes, "WWF_ecoregions"))) {
  download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
                destfile = file.path(Dir.Shapes, "wwf_ecoregions.zip"))
  unzip(file.path(Dir.Shapes, "wwf_ecoregions.zip"), exdir = file.path(Dir.Shapes, "WWF_ecoregions"))
}

EcoregionMask <- read_sf(file.path(Dir.Base,"data","shapes", "WWF_ecoregions", "official", "wwf_terr_ecos.shp"))
EcoregionMask <- st_make_valid(EcoregionMask)

EcoregionMask_hex <- st_make_valid(EcoregionMask)%>%
  st_transform("+proj=eqearth +wktext")

d1 <- read_csv("./data/organized_data_to_append/global_lake_res_CH4_emission_DB_all.csv") %>%
  #filter(data_source != "Malerba et al 2022") %>%
  group_by(lat, lon, ref, obs_year, country, continent, num_months_sampled, waterbody_type) %>%
  summarize(mean_ebu = mean(ch4_ebu, na.rm = T),
            sd_ebu = sd(ch4_ebu, na.rm = T),
            mean_diff = mean(ch4_diff, na.rm = T),
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

k <- d1 %>% ungroup(.) %>%
  select(lon, lat)

g <- d1 %>% ungroup(.) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")

link <- g %>%
  cbind(EcoregionMask_hex[st_nearest_feature(g, EcoregionMask_hex),]) %>%
  mutate(dist = st_distance(geometry, geometry.1, by_element = T)) %>%
  mutate(biome_type = case_when(
    BIOME == 1 ~ "TROPICAL MOIST FOREST",
    BIOME == 2 ~ "TROPICAL DRY FOREST",
    BIOME == 3 ~ "TROPICAL CONIFEROUS FOREST",
    BIOME == 4 ~ "TEMPERATE BROADLEAF FOREST",
    BIOME == 5 ~ "TEMPERATE CONIFEROUS FOREST",
    BIOME == 6 ~ "BOREAL FOREST",
    BIOME == 7 ~ "TROPICAL GRASSLAND",
    BIOME == 8 ~ "TEMPERATE GRASSLAND",
    BIOME == 9 ~ "FLOODED GRASSLAND",
    BIOME == 10 ~ "MONTANE GRASSLAND",
    BIOME == 11 ~ "TUNDRA",
    BIOME == 12 ~ "MEDITERRANIAN FOREST",
    BIOME == 13 ~ "DESERT",
    BIOME == 14 ~ "MANGROVES",
    BIOME == 98 ~ "LAKE",
    BIOME == 99 ~ "ROCK & ICE",
    TRUE ~ NA_character_)) %>%
  select(-OBJECTID, -AREA, -PERIMETER, -ECO_NAME, -REALM, -geometry.1, -dist, -eco_code, -PER_area, -PER_area_1,
         -PER_area_2, -area_km2, -Shape_Area, -Shape_Leng, -G200_STAT, -GBL_STAT, -G200_REGIO, -G200_NUM, -G200_BIOME,
         -ECO_SYM, -ECO_ID, -ECO_NUM, -BIOME) %>%
  bind_cols(., k) %>%
  write_csv(., "./data/organized_data_to_append/global_lake_res_DB_combined_refs_WWF_add.csv")
