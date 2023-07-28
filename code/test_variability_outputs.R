library(tidyverse)

jp <- ggplot2::map_data('world2', 'japan')
class(jp)


time_series <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/J_ebullition_variability.csv") %>%
  select(date, FOA_mean_ebullition, FOA_time_low_ebu, FOA_time_high_ebu, FOA_space_low_ebu, FOA_space_high_ebu,
         FOA_coefficient_low_ebu, FOA_coefficient_high_ebu) %>%
  reshape2::melt(., id="date") %>%
  group_by(date, variable) %>%
  summarise(mean.x = quantile(value, probs = 0.50),
            lower.x = quantile(value, probs = 0.05),
            upper.x = quantile(value, probs = 0.95))

ggplot(time_series, aes(date, mean.x, color = variable))+
  geom_line()

map <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/J_ebullition_variability.csv") %>%
  select(lat, lon, hylak_id, FOA_mean_ebullition, FOA_time_low_ebu, FOA_time_high_ebu, FOA_space_low_ebu, FOA_space_high_ebu,
         FOA_coefficient_low_ebu, FOA_coefficient_high_ebu, FOA_mean_error_low_ebu, FOA_mean_error_high_ebu) %>%
  reshape2::melt(., id=c("lat","lon","hylak_id")) %>%
  group_by(lat, lon,hylak_id, variable) %>%
  summarise(mean.x = mean(value))

map <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/J_ebullition_variability.csv") %>%
  select(lat, lon, hylak_id, FOA_mean_ebullition, FOA_time_low_ebu, FOA_time_high_ebu, FOA_space_low_ebu, FOA_space_high_ebu,
         FOA_coefficient_low_ebu, FOA_coefficient_high_ebu, FOA_mean_error_low_ebu, FOA_mean_error_high_ebu) %>%
  reshape2::melt(., id=c("lat","lon","hylak_id")) %>%
  group_by(lat, lon,hylak_id, variable) %>%
  summarise(mean.x = quantile(value, probs = 0.50))


library(ggfortify)
jp <-  map('world2', 'japan', plot = FALSE, fill = TRUE)
class(jp)


p <- autoplot(jp, geom = 'polygon', fill = "white")

ggplot(map) + 
  geom_point(aes(lon, lat, color = log(mean.x+1)), inherit.aes = F, pch = 16, 
               size = 1)+
  scale_color_gradient(low="darkblue",
                       high="red", space ="Lab", na.value="grey",
                       name = "**log(mg CH4 m-2 d-1)** ",
                       limits = c(-1,10)) +
  facet_wrap(~variable)

  
map_sums <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/J_ebullition_variability.csv") %>%
  select(lat, lon, hylak_id, FOA_mean_ebullition, FOA_time_low_ebu, FOA_time_high_ebu, FOA_space_low_ebu, FOA_space_high_ebu,
         FOA_coefficient_low_ebu, FOA_coefficient_high_ebu, FOA_mean_error_low_ebu, FOA_mean_error_high_ebu) %>%
  reshape2::melt(., id=c("lat","lon","hylak_id")) %>%
  group_by(variable) %>%
  summarise(mean.x = mean(value))
