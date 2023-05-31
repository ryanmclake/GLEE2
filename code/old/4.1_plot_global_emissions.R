library(vroom, warn.conflicts = FALSE)
library(sf, warn.conflicts = FALSE)
library(ggplot2, warn.conflicts = FALSE)
library(maps, warn.conflicts = FALSE)
library(hexbin, warn.conflicts = FALSE)
library(rnaturalearth, warn.conflicts = FALSE)

output_emission <- vroom::vroom("./output/global_ebullition_emission.csv", col_names = T) %>%
  na.omit(.) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") 

output_emission_year

output_emission_month

output_emission_all

world <-  ne_download(scale = 110, type = 'land', category = 'physical', returnclass = "sf") %>%
  st_transform("+proj=eqearth +wktext")

grid_spacing <- 100000 # CRS units in meters (100000 m = 111 km & 111 km ~ 1 Decimal degree)

years <- c(unique(output_emission$year))
months <- c(unique(output_emission$month))

for(i in 1:length(years)){
  for(g in 1:length(months)){
    
    output_emission <- output_emission %>%
      filter(year  == years[i]) %>%
      filter(month == months[g])
    
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
    
    area_hexes <- st_join(output_emission, grid, join = st_intersects)
    
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