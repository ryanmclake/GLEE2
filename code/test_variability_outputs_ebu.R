library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(patchwork)

world <-  ne_download(scale = 110, type = 'land', category = 'physical', returnclass = "sf") %>%
  st_transform("+proj=eqearth +wktext")

map_base <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_ebullition_variability_BASE_sum.csv") %>%
  mutate(`Baseline Ebullition Predictions` = ifelse(is.infinite(`Baseline Ebullition Predictions`), NA, `Baseline Ebullition Predictions`)) %>%
  mutate(`Baseline Ebullition Predictions` =`Baseline Ebullition Predictions`*0.001*365*0.4) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Baseline Ebullition Predictions` <= quantile(.$`Baseline Ebullition Predictions`, 0.95)) %>%
  filter(`Baseline Ebullition Predictions` >= quantile(.$`Baseline Ebullition Predictions`, 0.05)) %>%
  #mutate(`Baseline Ebullition Predictions` = log(`Baseline Ebullition Predictions`*0.001*365*0.4+1)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Baseline Ebullition Predictions`))+
  geom_sf(data = world, lwd =0.5, color = "black", fill = NA)+
  #geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5)+
  labs(title = "  Baseline Calibration")+
  scale_color_gradient(low="blue",
                       high="red", space ="Lab", na.value="black",
                       name = "**Ebullition Flux** <br> log(g CH<sub>4</sub> m<sup>-2</sup> yr<sup>-1</sup>)", limits = c(0,8.5)) +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-7000000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = c(0.11, 0.35),
        legend.direction = "vertical",
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = element_text(size=13),
        legend.key.height  = unit(.5, 'cm'),
        legend.key.width =  unit(.3, 'cm'))


map_time_low <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_ebullition_variability_TIME_sum.csv") %>%
  mutate(`Time-low Ebu Predictions` = ifelse(is.infinite(`Time-low Ebu Predictions`), NA, `Time-low Ebu Predictions`)) %>%
  mutate(`Time-low Ebu Predictions` =`Time-low Ebu Predictions`*0.001*365*0.4) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Time-low Ebu Predictions` <= quantile(.$`Time-low Ebu Predictions`, 0.95)) %>%
  filter(`Time-low Ebu Predictions` >= quantile(.$`Time-low Ebu Predictions`, 0.05)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Time-low Ebu Predictions`))+
  geom_sf(data = world, lwd =0.5, color = "black", fill = NA)+
  #geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5)+
  labs(title = "  Temporal Error (-1 S.D.)")+
  scale_color_gradient(low="blue",
                       high="red", space ="Lab", na.value="black",
                       name = "**Measurement Error (Time-low)** <br> log(g CH<sub>4</sub> m<sup>-2</sup> yr<sup>-1</sup>)", limits = c(0,8.5)) +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-7000000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = "none")


map_time_high <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_ebullition_variability_TIME_sum.csv") %>%
  mutate(`Time-high Ebu Predictions` = ifelse(is.infinite(`Time-high Ebu Predictions`), NA, `Time-high Ebu Predictions`)) %>%
  mutate(`Time-high Ebu Predictions` =`Time-high Ebu Predictions`*0.001*365*0.4) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Time-high Ebu Predictions` <= quantile(.$`Time-high Ebu Predictions`, 0.95)) %>%
  filter(`Time-high Ebu Predictions` >= quantile(.$`Time-high Ebu Predictions`, 0.05)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Time-high Ebu Predictions`))+
  geom_sf(data = world, lwd =0.5, color = "black", fill = NA)+
  #geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5)+
  labs(title = "  Temporal Error (+1 S.D.)")+
  scale_color_gradient(low="blue",
                       high="red", space ="Lab", na.value="black",
                       name = "**Measurement Error (Time-high)** <br>log(mg CH4 m-1 yr-1)", limits = c(0,8.5)) +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-7000000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = "none")


map_space_low <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_ebullition_variability_SPACE_sum.csv") %>%
  mutate(`Space-low Ebu Predictions` = ifelse(is.infinite(`Space-low Ebu Predictions`), NA, `Space-low Ebu Predictions`)) %>%
  mutate(`Space-low Ebu Predictions` =`Space-low Ebu Predictions`*0.001*365*0.4) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Space-low Ebu Predictions` <= quantile(.$`Space-low Ebu Predictions`, 0.95)) %>%
  filter(`Space-low Ebu Predictions` >= quantile(.$`Space-low Ebu Predictions`, 0.05)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Space-low Ebu Predictions`))+
  geom_sf(data = world, lwd =0.5, color = "black", fill = NA)+
  #geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5)+
  labs(title = "  Spatial Error (-1 S.D.)")+
  scale_color_gradient(low="blue",
                       high="red", space ="Lab", na.value="black",
                       name = "**Measurement Error (Space-low)** <br>log(mg CH4 m-1 yr-1)", limits = c(0,8.5)) +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-7000000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = "none")


map_space_high <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_ebullition_variability_SPACE_sum.csv") %>%
  mutate(`Space-high Ebu Predictions` = ifelse(is.infinite(`Space-high Ebu Predictions`), NA, `Space-high Ebu Predictions`)) %>%
  mutate(`Space-high Ebu Predictions` =`Space-high Ebu Predictions`*0.001*365*0.4) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Space-high Ebu Predictions` <= quantile(.$`Space-high Ebu Predictions`, 0.95)) %>%
  filter(`Space-high Ebu Predictions` >= quantile(.$`Space-high Ebu Predictions`, 0.05)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Space-high Ebu Predictions`))+
  geom_sf(data = world, lwd =0.5, color = "black", fill = NA)+
  #geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5)+
  labs(title = "  Measurement Error-Space (+1 S.D.)")+
  scale_color_gradient(low="blue",
                       high="red", space ="Lab", na.value="black",
                       name = "**Measurement Error (Space-high)** <br>log(mg CH4 m-1 yr-1)", limits = c(0,8.5)) +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-7000000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = "none")



map_param_low <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_ebullition_variability_PARAM_sum.csv") %>%
  mutate(`Param-low Ebu Predictions` = ifelse(is.infinite(`Param-low Ebu Predictions`), NA, `Param-low Ebu Predictions`)) %>%
  mutate(`Param-low Ebu Predictions` =`Param-low Ebu Predictions`*0.001*365*0.4) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Param-low Ebu Predictions` <= quantile(.$`Param-low Ebu Predictions`, 0.95)) %>%
  filter(`Param-low Ebu Predictions` >= quantile(.$`Param-low Ebu Predictions`, 0.05)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Param-low Ebu Predictions`))+
  geom_sf(data = world, lwd =0.5, color = "black", fill = NA)+
  #geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5)+
  labs(title = "  Coefficient Error (-1 S.D.)")+
  scale_color_gradient(low="blue",
                       high="red", space ="Lab", na.value="black",
                       name = "**Parameter Error (low)** <br>log(mg CH4 m-1 yr-1)", limits = c(0,8.5)) +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-7000000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = "none")

map_param_high <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_ebullition_variability_PARAM_sum.csv") %>%
  mutate(`Param-high Ebu Predictions` = ifelse(is.infinite(`Param-high Ebu Predictions`), NA, `Param-high Ebu Predictions`)) %>%
  mutate(`Param-high Ebu Predictions` =`Param-high Ebu Predictions`*0.001*365*0.4) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Param-high Ebu Predictions` <= quantile(.$`Param-high Ebu Predictions`, 0.95)) %>%
  filter(`Param-high Ebu Predictions` >= quantile(.$`Param-high Ebu Predictions`, 0.05)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Param-high Ebu Predictions`))+
  geom_sf(data = world, lwd =0.5, color = "black", fill = NA)+
  #geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5)+
  labs(title = "  Coefficient Error (+1 S.D.)")+
  scale_color_gradient(low="blue",
                       high="red", space ="Lab", na.value="black",
                       name = "**Parameter Error (high)** <br>log(mg CH4 m-1 yr-1)", limits = c(0,8.5)) +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-7000000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = "none")


map_model_low <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_ebullition_variability_MODEL_sum.csv") %>%
  mutate(`Model-low Ebu Predictions` = ifelse(is.infinite(`Model-low Ebu Predictions`), NA, `Model-low Ebu Predictions`)) %>%
  mutate(`Model-low Ebu Predictions` =`Model-low Ebu Predictions`*0.001*365*0.4) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Model-low Ebu Predictions` <= quantile(.$`Model-low Ebu Predictions`, 0.95)) %>%
  filter(`Model-low Ebu Predictions` >= quantile(.$`Model-low Ebu Predictions`, 0.05)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Model-low Ebu Predictions`))+
  geom_sf(data = world, lwd =0.5, color = "black", fill = NA)+
  #geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5)+
  labs(title = "  Model Error (-1 S.D.)")+
  scale_color_gradient(low="blue",
                       high="red", space ="Lab", na.value="black",
                       name = "**Model Error (low)** <br>log(mg CH4 m-1 yr-1)", limits = c(0,8.5)) +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-7000000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = "none")

map_model_high <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_ebullition_variability_MODEL_sum.csv") %>%
  mutate(`Model-high Ebu Predictions` = ifelse(is.infinite(`Model-high Ebu Predictions`), NA, `Model-high Ebu Predictions`)) %>%
  mutate(`Model-high Ebu Predictions` =`Model-high Ebu Predictions`*0.001*365*0.4) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Model-high Ebu Predictions` <= quantile(.$`Model-high Ebu Predictions`, 0.95)) %>%
  filter(`Model-high Ebu Predictions` >= quantile(.$`Model-high Ebu Predictions`, 0.05)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Model-high Ebu Predictions`))+
  geom_sf(data = world, lwd =0.5, color = "black", fill = NA)+
  #geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5)+
  labs(title = "  Model Error (+1 S.D.)")+
  scale_color_gradient(low="blue",
                       high="red", space ="Lab", na.value="black",
                       name = "**Model Error (high)** <br> mg CH4 m-1 yr-1", limits = c(0,8.5)) +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-7000000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = "none")

maps <- (map_base + map_time_low + map_time_high)/
        (map_space_low + map_space_high + map_model_low)/
        (map_model_high + map_param_low + map_param_high)



maps









line <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_ebullition_variability_BASE_mean.csv") %>%
  mutate(`Baseline Ebullition Predictions` = log(`Baseline Ebullition Predictions`+1)) %>%
  mutate(centr_lat = ceiling(centr_lat)) %>%
  group_by(centr_lat) %>%
  summarise(`Baseline Ebullition Predictions` = mean(`Baseline Ebullition Predictions`))

plot(line$`Baseline Ebullition Predictions`, line$centr_lat, type = "l", col = "red")
  



# Set the grid sizing to overlay on the world data

# CRS units in meters (100000 m = 111 km & 111 km ~ 1 Decimal degree) ??? Not sure of this!!!
cs <- c(100000, 100000)

# Set up our spatial grid
grid <- st_make_grid(
  world,
  cellsize = cs,
  #n = c(200, 200), # grid granularity
  crs = st_crs(world),
  what = "polygons",
  flat_topped = T,
  square = F) %>%
  st_intersection(world)


# assign an index column in our grid that will match up with our area hexes (below)
grid <- st_sf(index = 1:length(lengths(grid)), grid)

# Join all of the global slope data with the grid overlay
area_hexes <- st_join(map_base, grid, join = st_intersects)

# Calculate the median slope among all of the lakes that are in the respective grid bins
area_hexes_avg <- area_hexes %>%
  st_drop_geometry() %>%
  group_by(index) %>%
  mutate(bin_count = n()) %>%
  summarise(`Baseline Ebullition Predictions` = mean(`Baseline Ebullition Predictions`, na.rm = TRUE),
            bin_count = median(bin_count)) %>%
  right_join(grid, by="index") %>%
  st_sf()

# make the global plot
flux_map <-
  ggplot() +
  geom_sf(data = map_base, lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Baseline Ebullition Predictions`))+
  geom_sf(data = world, lwd =0.5, color = "black", fill = NA)+
  #geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5)+
  scale_color_gradient(low="blue",
                       high="red", space ="Lab", na.value="black",
                       name = "**Ebullition Flux** <br>log(mg CH4 m-1 yr-1)") +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-8600000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = c(0.11, 0.35),
        legend.direction = "vertical",
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = element_text(size=9),
        legend.key.height  = unit(.5, 'cm'),
        legend.key.width =  unit(.3, 'cm'))

flux_map


# make the global plot
flux_map2 <-
  ggplot() +
  geom_sf(data = world, lwd = 0.5, color = "black")+
  geom_sf(data = area_hexes_avg,lwd = 0.05,
          aes(fill = `Baseline Ebullition Predictions`))+
  #geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5)+
  scale_fill_gradient2(midpoint=4.2, low="blue", mid="yellow",
                        high="red", space ="Lab", na.value="grey",
                        name = "**Interannual Ebu Rate** <br>(mg CH4 m-1 yr-1)") +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-8600000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = c(0.11, 0.35),
        legend.direction = "vertical",
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = element_text(size=9),
        legend.key.height  = unit(.5, 'cm'),
        legend.key.width =  unit(.3, 'cm'))

flux_map2


