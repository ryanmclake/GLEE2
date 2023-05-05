
library(dplyr)
library(vroom)
library(sf)
library(dplyr)
library(ggplot2)
library(maps)
library(patchwork)
library(hexbin)
library(rnaturalearth)
library(units)
library(zoo)
library(randomForest)
library(viridis)
library(rpart)
library(rpart.plot)
library(grid)
library(outliers)
library(trend, warn.conflicts = FALSE)
library(broom, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
library(Kendall)
library(tidyr)

set.seed(32)

FOA_model <- function(temp){
  rnorm(1, 96.98125, 17.46113) * rnorm(1, 1.17077, 0.02826) ^ (temp-20)
}

Aben_model <- function(temp){
  rnorm(1, 100, 1) * rnorm(1, 1.1, 0.01) ^ (temp-20)
}

dat <- vroom("./data/Germany.csv", col_names = F) %>%
  filter(X1 >= 2000)

out <- list()

for(i in 1:10){
  
  dat2 <- dat %>% 
    na.omit(.) %>%
    mutate(ebullition_first_order = ifelse(X9 > 50, 0, mapply(FOA_model, X8-273.15)),
           ebullition_aben_eqn = ifelse(X9 > 50, 0, mapply(Aben_model, X8-273.15)),
           ebullition_Lit_area_corrected_first_order = ebullition_first_order*X7*1000000*0.000001*0.4*30,
           ebullition_Lit_area_corrected_aben_eqn = ebullition_aben_eqn*X7*1000000*0.000001*0.4*30)
  
  out[[i]] <- dat2
  
}

out = do.call(rbind, out) 

out_sf <- out %>%
  st_as_sf(coords = c("X5", "X4"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") 

out2 <- out %>%  group_by(X1, X2, X3, X4, X5) %>%
  summarize(ebullition_first_order_mean = mean(ebullition_first_order),
            ebullition_first_order_sd = sd(ebullition_first_order),
            ebullition_aben_eqn_mean = mean(ebullition_aben_eqn),
            ebullition_aben_eqn_sd = sd(ebullition_aben_eqn),
            ebullition_Lit_area_corrected_first_order_mean = mean(ebullition_Lit_area_corrected_first_order),
            ebullition_Lit_area_corrected_first_order_sd = sd(ebullition_Lit_area_corrected_first_order),
            ebullition_Lit_area_corrected_aben_eqn_mean = mean(ebullition_Lit_area_corrected_aben_eqn),
            ebullition_Lit_area_corrected_aben_eqn_sd = sd(ebullition_Lit_area_corrected_aben_eqn)) %>%
  st_as_sf(coords = c("X5", "X4"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") 

grid_spacing <- 10000 # CRS units in meters (100000 m = 111 km & 111 km ~ 1 Decimal degree)

grid <- st_make_grid(
  world,
  cellsize = c(grid_spacing, grid_spacing),
  #n = c(200, 200), # grid granularity
  crs = st_crs(world),
  what = "polygons",
  flat_topped = T,
  square = F) %>%
  st_intersection(world)

grid <- st_sf(index = 1:length(lengths(grid)), grid)

area_hexes <- st_join(out_sf, grid, join = st_intersects)
  
area_hexes_avg <- area_hexes %>%
  st_drop_geometry() %>%
  group_by(index, X2) %>%
  mutate(bin_count = n()) %>%
  summarize(ebullition_first_order_mean = mean(ebullition_first_order),
            ebullition_first_order_sd = sd(ebullition_first_order),
            ebullition_aben_eqn_bin = mean(ebullition_aben_eqn),
            ebullition_aben_eqn_sd = sd(ebullition_aben_eqn),
            ebullition_Lit_area_corrected_first_order_mean = mean(ebullition_Lit_area_corrected_first_order),
            ebullition_Lit_area_corrected_first_order_sd = sd(ebullition_Lit_area_corrected_first_order),
            ebullition_Lit_area_corrected_aben_eqn_mean = mean(ebullition_Lit_area_corrected_aben_eqn),
            ebullition_Lit_area_corrected_aben_eqn_sd = sd(ebullition_Lit_area_corrected_aben_eqn),
            bin_count = median(bin_count)) %>%
  right_join(grid, by="index") %>%
  st_sf()

ebu_base_temp_compare <- base %>% 
  select(ch4_ebu, country, month, year, tot_sampling_events) %>%
  filter(country == "Germany") %>%
  na.omit(.) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events <= 10, 29.434313, NA)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 10, 8.903169, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 20, 1.914520, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 30, 1.131719, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 40, 1.123484, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 50, 1.101714, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 60, 1.154789, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 100, 1.151526, ebu_sd))

month_ts <- area_hexes_avg %>% 
  ggplot(.)+
  geom_line(aes(x = X2, y = ebullition_first_order_mean, group = index), lwd = 0.1, alpha = 0.1, color = "blue")+
  geom_ribbon(aes(x = X2, y = ebullition_first_order_mean,
                  ymin = ebullition_first_order_mean-ebullition_first_order_sd,
                  ymax = ebullition_first_order_mean+ebullition_first_order_sd,
                  group = index),
              alpha = 0.005, fill = "midnightblue") +
  geom_point(data = ebu_base_temp_compare, aes(x = month, y = ch4_ebu), color = "red", size = 2)+
  geom_errorbar(data = ebu_base_temp_compare, aes(x = month, y = ch4_ebu, ymin = ch4_ebu-ebu_sd, ymax = ch4_ebu+ebu_sd),width = 0.2)+
  ylab("Ebullition Rate (mg CH4 m-2 d-1)")+
  xlab("Month")+
  theme_classic()

world <-  ne_download(scale = 110,  type = 'countries', returnclass = "sf") %>%
  st_transform("+proj=eqearth +wktext") %>%
  filter(SOVEREIGNT == "Germany")

grid_spacing <- 50000 # CRS units in meters (100000 m = 111 km & 111 km ~ 1 Decimal degree)

years <- c(unique(dat2$X1))
months <- c(unique(dat2$X2))

for(i in 1:length(years)){
  for(g in 1:length(months)){
    
    dat3 <- dat2 %>%
      filter(X1  == years[i]) %>%
      filter(X2 == months[g])%>%
      st_as_sf(coords = c("X5", "X4"), crs = 4326) %>%
      st_transform("+proj=eqearth +wktext")
    
    grid <- st_make_grid(
      world,
      cellsize = c(grid_spacing, grid_spacing),
      #n = c(200, 200), # grid granularity
      crs = st_crs(world),
      what = "polygons",
      flat_topped = T,
      square = F) %>%
      st_intersection(world)
    
    grid <- st_sf(index = 1:length(lengths(grid)), grid)
    
    area_hexes <- st_join(dat3, grid, join = st_intersects)
    
    area_hexes_avg <- area_hexes %>%
      st_drop_geometry() %>%
      group_by(X1, X2, index) %>%
      mutate(bin_count = n()) %>%
      summarise(bin_emissions_FO = exp(mean(log((ebullition_first_order)))),
                bin_lake_littoral_area_m2 = sum(X16*1000000*0.3),
                bin_count = median(bin_count)) %>%
      mutate(littoral_area_weighted_emission = bin_emissions_FO*bin_lake_littoral_area_m2*0.000001) %>%
      right_join(grid, by="index") %>%
      st_sf()
    
    p <- ggplot() +
      geom_sf(data = world, lwd = 0.5, color = "black")+
      geom_sf(data = area_hexes_avg,lwd = 0.05,
              aes(fill = littoral_area_weighted_emission))+
      labs(title = paste0("Spains's littoral area weighted monthly ebullition rate (year:",years[i],"  month:", months[g],")"))+
      geom_sf_text(data = area_hexes_avg, aes(label = bin_count), size = 1.5, color = "white")+
      scale_fill_gradient(low="black", high="turquoise1", space ="Lab", na.value="grey",
                          name = "Ebullition (kg CH4/month)", limits=c(0,100000)) +
      #coord_sf(xlim = c(-15000000, 16000000), ylim = c(-8600000, 8600000), expand = FALSE) +
      guides(fill = guide_colourbar(title.position = "top"))+
      theme_void()+
      theme(legend.position = c(0.05, 0.25),
            legend.direction = "vertical",
            legend.title = ggtext::element_markdown(size = 10),
            legend.text = element_text(size=9),
            legend.key.height  = unit(.5, 'cm'),
            legend.key.width =  unit(.3, 'cm'))
    
    ggsave(p, path = ".",filename = paste0("./spain_",years[i],"_",months[g],".jpg"),
                  width = 8, height = 4, device='jpg', dpi=175)
    
  }
}

  
