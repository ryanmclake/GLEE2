
boxplotit<-function(data){
  
  n <- boxplot(log(value+1) ~ variable, data = data, 
          outline = F, 
          ylab = "log(mg CH4 m-2 d-1)",
          main = paste0("Ebullition Flux Partition for ", months[g]), 
          ylim = c(0,15))
  ratio <- tibble::tibble(a = c(diff(n$stats[c(2,4), 1])*0.1289776,diff(n$stats[c(2,4), 1])*0.8710224))
  rect(xleft = c(0.6) + seq_along(n$n[1])-1, xright = 1.4 + seq_along(n$n[1])-1, ybottom = n$stats[2, 1]+ratio$a[1]+ratio$a[2],ytop = n$stats[2, 1], col = 2)
  rect(xleft = c(0.6) + seq_along(n$n[1])-1, xright = 1.4 + seq_along(n$n[1])-1, ybottom = n$stats[2, 1], ytop = n$stats[2, 1]+ratio$a[1], col = 1)
  legend("topleft",                    # Add legend to plot
         legend = c("Model Uncertainty", "Parameter Uncertainty"),
         # col = c("#56B4E9":"#D55E00"),
         col = 2:1,
         pch = 15,
         cex = 0.8)
}




lake <- vroom::vroom("./output/global_lake_emissions 1.csv", delim = " ")

ggplot(lake, aes(ebu_Aben_estimate, mean))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1, color = "red", lwd = 2, alpha = 0.5)+
  labs(title = paste0("Lake Emission Comparisons (N = ",length(unique(lake$hylak_id)),")"))





res <- vroom::vroom("./output/global_reservoir_emissions 1.csv", delim = " ") %>%
  na.omit(.) %>%
  mutate(outlier = (abs(mean - median(mean)) > 2*sd(mean)) %>% as.vector()) %>%
  # turn those outliers into an NA value
  mutate(mean = ifelse(outlier == "TRUE",NA,mean))%>%
  na.omit(.)
  mutate(outlier = (abs(mean - median(mean)) > 2*sd(mean)) %>% as.vector()) %>%
  # turn those outliers into an NA value
  mutate(mean = ifelse(outlier == "TRUE",NA,mean))%>%
  na.omit(.)


res_avg <- res %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  mutate(month = lubridate::month(date))


world <-  ne_download(scale = 110, type = 'land', category = 'physical', returnclass = "sf") %>%
  st_transform("+proj=eqearth +wktext")

# Set the grid sizing to overlay on the world data
grid_spacing <- 300000 # CRS units in meters (100000 m = 111 km & 111 km ~ 1 Decimal degree) ??? Not sure of this!!!

# Set up our spatial grid
grid <- st_make_grid(
  world,
  cellsize = c(grid_spacing, grid_spacing),
  #n = c(200, 200), # grid granularity
  crs = st_crs(world),
  what = "polygons",
  flat_topped = T,
  square = F) %>%
  st_intersection(world)

# assign an index column in our grid that will match up with our area hexes (below)
grid <- st_sf(index = 1:length(lengths(grid)), grid)

# Join all of the global slope data with the grid overlay
area_hexes <- st_join(res_avg, grid, join = st_intersects)

# Calculate the median slope among all of the lakes that are in the respective grid bins
area_hexes_avg_mean <- area_hexes %>%
  ungroup(.) %>%
  st_drop_geometry() %>%
  group_by(date, index) %>%
  mutate(bin_count = n()) %>%
  summarise(variance = median(abs(sd_low - sd)/bin_count, na.rm = TRUE),
            mean_flux = median(mean, na.rm = TRUE),
            mean_flux_low = median((mean-sd_low), na.rm = TRUE),
            mean_flux_high = median(mean+sd, na.rm = TRUE),
            Aben_mean_flux = median(ebu_Aben_estimate, na.rm = TRUE),
            Aben_sd_flux = sd(ebu_Aben_estimate, na.rm = TRUE),
            bin_count = median(bin_count)) %>%
  right_join(grid, by="index") %>%
  st_sf()

summary(area_hexes_avg_mean$variance)

months <- c(unique(area_hexes_avg_mean$date))

for(g in 1: length(months)){
# make the global plot
plot <- area_hexes_avg_mean %>%
  filter(date == months[g]) %>%
  ggplot(.) +
  geom_sf(data = world, lwd = 0.5, color = "black")+
  geom_sf(lwd = 0.05,
          aes(fill = log(variance+1)))+
  #geom_sf_text(data = area_hexes_avg_mean, aes(label = bin_count), size = 1)+
  labs(title = paste0(".   ",months[g]))+
  scale_fill_gradient(low="grey",
                       high="blue4", space ="Lab", na.value="white",
                       name = "**Monthly bin error** <br> mg CH4 m-2 d-1",
                      limits=c(0,8)) +
  coord_sf(xlim = c(-15000000, 16000000), ylim = c(-8600000, 8600000), expand = FALSE) +
  guides(fill = guide_colourbar(title.position = "top"))+
  theme_void()+
  theme(legend.position = c(0.11, 0.35),
        legend.direction = "vertical",
        legend.title = ggtext::element_markdown(size = 10),
        legend.text = element_text(size=9),
        legend.key.height  = unit(.5, 'cm'),
        legend.key.width =  unit(.3, 'cm'))

boxplot <- area_hexes_avg_mean %>%
  st_drop_geometry() %>%
  select(-mean_flux_low, -mean_flux_high, -bin_count, -index) %>%
  filter(date == months[g])%>%
  reshape2::melt(., id.vars = c("date"))

jpeg(paste0("./output/figures/month_",months[g],"_global_CH4_ebullition_boxplot.jpg"), width = 300, height = 300)
boxplotit(data = boxplot)
dev.off()

# Save it as a JPEG picture file (THIS IS A REALLY HIGH RESOLUTION) --> almost 50 MB figure
ggsave(plot, path = ".",
       filename = paste0("./output/figures/month_",months[g],"_global_CH4_ebullition.jpg"),
       width = 8, height = 6, device='jpg', dpi=300)

}




# ggplot(res, aes(ebu_Aben_estimate, mean, color = mean_temp_k), alpha = 0.5)+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1, color = "red", lwd = 2, alpha = 0.5)+
#   labs(title = paste0("Lake Emission Comparisons (N = ",length(unique(lake$hylak_id)),")"))


