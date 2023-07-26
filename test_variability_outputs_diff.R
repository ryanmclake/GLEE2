library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(patchwork)

world <-  ne_download(scale = 110, type = 'land', category = 'physical', returnclass = "sf") %>%
  st_transform("+proj=eqearth +wktext")

map_base <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_diffusion_variability_BASE_sum.csv") %>%
  mutate(`Baseline Diffusion Predictions` = ifelse(is.infinite(`Baseline Diffusion Predictions`), NA, `Baseline Diffusion Predictions`)) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Baseline Diffusion Predictions` <= quantile(.$`Baseline Diffusion Predictions`, 0.95)) %>%
  filter(`Baseline Diffusion Predictions` >= quantile(.$`Baseline Diffusion Predictions`, 0.05)) %>%
  mutate(`Baseline Diffusion Predictions` = log(`Baseline Diffusion Predictions`*0.001*365+1)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Baseline Diffusion Predictions`))+
  geom_sf(data = world, lwd =0.5, color = "black", fill = NA)+
  #geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5)+
  labs(title = "  Baseline Calibration")+
  scale_color_gradient(low="blue",
                       high="red", space ="Lab", na.value="black",
                       name = "**Diffusion Flux** <br> log(g CH<sub>4</sub> m<sup>-2</sup> yr<sup>-1</sup>)", limits = c(0,8.5)) +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-7000000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = c(0.11, 0.35),
        legend.direction = "vertical",
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = element_text(size=13),
        legend.key.height  = unit(.5, 'cm'),
        legend.key.width =  unit(.3, 'cm'))


map_time_low <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_diffusion_variability_TIME_sum.csv") %>%
  mutate(`Time-low Diff Predictions` = ifelse(is.infinite(`Time-low Diff Predictions`), NA, `Time-low Diff Predictions`)) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Time-low Diff Predictions` <= quantile(.$`Time-low Diff Predictions`, 0.95)) %>%
  filter(`Time-low Diff Predictions` >= quantile(.$`Time-low Diff Predictions`, 0.05)) %>%
  mutate(`Time-low Diff Predictions` = log(`Time-low Diff Predictions`*0.001*365+1)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Time-low Diff Predictions`))+
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


map_time_high <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_diffusion_variability_TIME_sum.csv") %>%
  mutate(`Time-high Diff Predictions` = ifelse(is.infinite(`Time-high Diff Predictions`), NA, `Time-high Diff Predictions`)) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Time-high Diff Predictions` <= quantile(.$`Time-high Diff Predictions`, 0.95)) %>%
  filter(`Time-high Diff Predictions` >= quantile(.$`Time-high Diff Predictions`, 0.05)) %>%
  mutate(`Time-high Diff Predictions` = log(`Time-high Diff Predictions`*0.001*365+1)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Time-high Diff Predictions`))+
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


map_space_low <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_diffusion_variability_SPACE_sum.csv") %>%
  mutate(`Space-low Diff Predictions` = ifelse(is.infinite(`Space-low Diff Predictions`), NA, `Space-low Diff Predictions`)) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Space-low Diff Predictions` <= quantile(.$`Space-low Diff Predictions`, 0.95)) %>%
  filter(`Space-low Diff Predictions` >= quantile(.$`Space-low Diff Predictions`, 0.05)) %>%
  mutate(`Space-low Diff Predictions` = log(`Space-low Diff Predictions`*0.001*365+1)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Space-low Diff Predictions`))+
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


map_space_high <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_diffusion_variability_SPACE_sum.csv") %>%
  mutate(`Space-high Diff Predictions` = ifelse(is.infinite(`Space-high Diff Predictions`), NA, `Space-high Diff Predictions`)) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Space-high Diff Predictions` <= quantile(.$`Space-high Diff Predictions`, 0.95)) %>%
  filter(`Space-high Diff Predictions` >= quantile(.$`Space-high Diff Predictions`, 0.05)) %>%
  mutate(`Space-high Diff Predictions` = log(`Space-high Diff Predictions`*0.001*365+1)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Space-high Diff Predictions`))+
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



map_param_low <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_diffusion_variability_PARAM_sum.csv") %>%
  mutate(`Param-low Diff Predictions` = ifelse(is.infinite(`Param-low Diff Predictions`), NA, `Param-low Diff Predictions`)) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Param-low Diff Predictions` <= quantile(.$`Param-low Diff Predictions`, 0.95)) %>%
  filter(`Param-low Diff Predictions` >= quantile(.$`Param-low Diff Predictions`, 0.05)) %>%
  mutate(`Param-low Diff Predictions` = log(`Param-low Diff Predictions`*0.001*365+1)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Param-low Diff Predictions`))+
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

map_param_high <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_diffusion_variability_PARAM_sum.csv") %>%
  mutate(`Param-high Diff Predictions` = ifelse(is.infinite(`Param-high Diff Predictions`), NA, `Param-high Diff Predictions`)) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Param-high Diff Predictions` <= quantile(.$`Param-high Diff Predictions`, 0.95)) %>%
  filter(`Param-high Diff Predictions` >= quantile(.$`Param-high Diff Predictions`, 0.05)) %>%
  mutate(`Param-high Diff Predictions` = log(`Param-high Diff Predictions`*0.001*365+1)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Param-high Diff Predictions`))+
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


map_model_low <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_diffusion_variability_MODEL_sum.csv") %>%
  mutate(`Model-low Diff Predictions` = ifelse(is.infinite(`Model-low Diff Predictions`), NA, `Model-low Diff Predictions`)) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Model-low Diff Predictions` <= quantile(.$`Model-low Diff Predictions`, 0.95)) %>%
  filter(`Model-low Diff Predictions` >= quantile(.$`Model-low Diff Predictions`, 0.05)) %>%
  mutate(`Model-low Diff Predictions` = log(`Model-low Diff Predictions`*0.001*365+1)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Model-low Diff Predictions`))+
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

map_model_high <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/Global_diffusion_variability_MODEL_sum.csv") %>%
  mutate(`Model-high Diff Predictions` = ifelse(is.infinite(`Model-high Diff Predictions`), NA, `Model-high Diff Predictions`)) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("centr_lon", "centr_lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(`Model-high Diff Predictions` <= quantile(.$`Model-high Diff Predictions`, 0.95)) %>%
  filter(`Model-high Diff Predictions` >= quantile(.$`Model-high Diff Predictions`, 0.05)) %>%
  mutate(`Model-high Diff Predictions` = log(`Model-high Diff Predictions`*0.001*365+1)) %>%
  ggplot(.) +
  geom_sf(lwd = 0.05, pch=16,size = 0.1,
          aes(color = `Model-high Diff Predictions`))+
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