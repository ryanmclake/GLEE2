#### Libraries #### 
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(vroom, warn.conflicts = FALSE)
library(sf, warn.conflicts = FALSE)
library(units, warn.conflicts = FALSE)
library(broom, warn.conflicts = FALSE)
library(Kendall, warn.conflicts = FALSE)
library(arrow, warn.conflicts = FALSE)
library(tidyverse, warn.conflicts = FALSE)
library(doParallel, warn.conflicts = FALSE)

s = Sys.time()

glee <- read_csv("./data/organized_data_to_append/global_lake_res_DB_combined_refs_WWF_sd_added.csv") %>% 
    select(-geometry) %>%
    st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
    st_transform("+proj=eqearth +wktext") %>%
    rename(year = obs_year)

lat_lon <- read_csv("./data/organized_data_to_append/global_lake_res_DB_combined_refs_WWF_sd_added.csv") %>% 
  select(lat, lon) 
  
glcp <- vroom::vroom("./data/GLCP/D1_glcp_slim_yearly_median2.csv") %>% 
  select(year,hybas_id,centr_lat,centr_lon,elevation,mean_annual_temp_k) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")

glcp_glee_join <- do.call('rbind', lapply(split(glee, 1:nrow(glee)), function(x) {
  st_join(x, glcp[glcp$year == unique(x$year),], join = st_nearest_feature)
})) %>% bind_cols(., lat_lon) %>% st_drop_geometry() %>%
  write.table(., file = paste0("./data/organized_data_to_append/global_lake_res_DB_refs_WWF_GLCP.txt"),
            append = T,
            row.names = T,
            col.names = !file.exists("./data/organized_data_to_append/global_lake_res_DB_refs_WWF_GLCP.txt"))

#### Time check ####
e <- Sys.time()
t=e-s
print(t)


