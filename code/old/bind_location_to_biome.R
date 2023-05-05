library(sf)
library(tidyverse)

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

k <- d3 %>% ungroup(.) %>%
  select(lon, lat)

d <- d3 %>% ungroup(.) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")

link <- d %>%
  cbind(EcoregionMask_hex[st_nearest_feature(d, EcoregionMask_hex),]) %>%
  mutate(dist = st_distance(geometry, geometry.1, by_element = T))

linked <- link %>%
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
  write_csv(., "./data/organized_data_to_append/global_lake_res_DB_with_WWF.csv")


