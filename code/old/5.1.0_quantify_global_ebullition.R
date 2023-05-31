# =======================================================================
#------------------------------------------------------------------------
set.seed(1236)

#### Libraries #### 
library(dplyr, warn.conflicts = FALSE)
library(units, warn.conflicts = FALSE)
library(doParallel, warn.conflicts = FALSE)
library(utils, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(vroom, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)

Aben_model_paper <- function(temp){
  est = 100 * 1.1 ^ (temp-20)
  return(est)
}

ebu_arhennius_lakes <- function(temp){
  est = rnorm(1,73.9, sd = 8.99) * rnorm(1,1.05, sd = 0.0209) ^ (temp-20) + rnorm(1,0, sd = 118)
  return(est)
}

ebu_arhennius_reservoirs <- function(temp){
  est = rnorm(1,38.4, sd = 32.3) * rnorm(1,1.54, sd = 0.185) ^ (temp-20) + rnorm(1,0, sd = 292) 
  return(est)
}

country <- c(list.files("/Users/ryanmcclure/Documents/GLEE2/data/countries"))

calc_global_flux <- function(country){

d_lake <- vroom::vroom(paste0("/Users/ryanmcclure/Documents/GLEE2/data/countries/",country), col_names = F, delim = ",") %>%
  rename(year = X1, month= X2,
         hylak_id = X3, lat = X4,
         lon = X5, continent = X6, country = X7,
         total_km2 = X8, waterbody_type = X9,
         mean_temp_k = X10, ice_cover_mean = X11) %>%
  mutate(ice_cover_mean = ifelse(is.na(ice_cover_mean),0,ice_cover_mean)) %>%
  filter(waterbody_type == 1) %>%
  na.omit(.)

#calculate_emission <- function(country){

d_aben_lake <- d_lake %>%
  group_by(hylak_id, lat, lon, total_km2) %>%
  mutate(ebu_Aben_estimate = ifelse(ice_cover_mean > 60, 0, mapply(Aben_model_paper, mean_temp_k-273.15)),
         'date' = lubridate::make_date(year = year, month = month)) %>%
  ungroup(.) %>%
  select(-year, -month)

out_jags <- list()

iteration <- seq(1:50)

for(i in 1:length(iteration)){
  
  d1 <- d_lake %>%
    group_by(hylak_id, lat, lon, total_km2) %>%
    mutate(ebu_JAGS_lake_estimate = ifelse(ice_cover_mean > 60, 0, mapply(ebu_arhennius_lakes, mean_temp_k-273.15)),
           'date' = lubridate::make_date(year = year, month = month),
           iteration = iteration[i])
  
  out_jags[[i]] <- d1
  
}

# e1 <- as.data.frame(do.call(rbind, out_jags))

# Raw Data to plot with shadows
e2 <- as.data.frame(do.call(rbind, out_jags)) %>%
  group_by(date, hylak_id) %>%
  summarize(mean = mean(ebu_JAGS_lake_estimate),
            sd = sd(ebu_JAGS_lake_estimate)) %>%
  mutate(mean=ifelse(mean<=0,0,mean)) %>%
  mutate(sd_low = ifelse(mean-sd < 0,0,sd),
         sd_low = ifelse(sd_low == 0, mean, sd_low)) %>%
  left_join(., d_aben_lake, by = c("date", "hylak_id")) %>%
  utils::write.table(., file = paste0("/Users/ryanmcclure/Documents/GLEE2/output/global_lake_emissions.csv"),
                     append = T,
                     row.names = F,
                     col.names = !file.exists("/Users/ryanmcclure/Documents/GLEE2/output/global_lake_emissions.csv"))

# p1 <- ggplot(e2, aes(x = date, y = mean))+
#   geom_line(lwd = 1, alpha = 1, color = "blue")+
#   geom_ribbon(aes(x = date, y = mean, ymin = mean-sd_low, ymax = mean+sd), alpha = 0.2, lwd = 0)+
#   geom_line(data = d_aben, aes(date, mean), lwd = 0.5, color = "red", inherit.aes = F, alpha = 0.3)+
#   coord_cartesian(ylim = c(-10,2000))+
#   labs(title = "Monthly Posterior Predictions of Japan's Lakes")+
#   theme_classic()

# The model calibration seems to be overestimating emissions from LAKES? 

# RESERVOIRS

d_res <- vroom::vroom(paste0("/Users/ryanmcclure/Documents/GLEE2/data/countries/",country), col_names = F, delim = ",") %>%
  rename(year = X1, month= X2,
         hylak_id = X3, lat = X4,
         lon = X5, continent = X6, country = X7,
         total_km2 = X8, waterbody_type = X9,
         mean_temp_k = X10, ice_cover_mean = X11) %>%
  mutate(ice_cover_mean = ifelse(is.na(ice_cover_mean),0,ice_cover_mean)) %>%
  filter(waterbody_type != 1) %>%
  na.omit(.)

#calculate_emission <- function(country){

d_aben_res <- d_res %>%
  group_by(hylak_id, lat, lon, total_km2, continent, country) %>%
  mutate(ebu_Aben_estimate = ifelse(ice_cover_mean > 60, 0, mapply(Aben_model_paper, mean_temp_k-273.15)),
         'date' = lubridate::make_date(year = year, month = month)) %>%
  ungroup(.) %>%
  select(-year, -month)

out_jags <- list()

iteration <- seq(1:50)

for(i in 1:length(iteration)){
  
  d1 <- d_res %>%
    group_by(hylak_id, lat, lon, total_km2) %>%
    mutate(ebu_JAGS_reservoir_estimate = ifelse(ice_cover_mean > 60, 0, mapply(ebu_arhennius_reservoirs, mean_temp_k-273.15)),
           'date' = lubridate::make_date(year = year, month = month),
           iteration = iteration[i])
  
  out_jags[[i]] <- d1
  
}

# e3 <- as.data.frame(do.call(rbind, out_jags))

# Raw Data to plot with shadows
e4 <- as.data.frame(do.call(rbind, out_jags)) %>%
  group_by(date, hylak_id) %>%
  summarize(mean = mean(ebu_JAGS_reservoir_estimate),
            sd = sd(ebu_JAGS_reservoir_estimate)) %>%
  mutate(mean=ifelse(mean<=0,0,mean)) %>%
  mutate(sd_low = ifelse(mean-sd < 0,0,sd),
         sd_low = ifelse(sd_low == 0, mean, sd_low),
         sd = ifelse(sd > 2000, mean, sd)) %>%
  left_join(., d_aben_res, by = c("date", "hylak_id")) %>%
  utils::write.table(., file = paste0("/Users/ryanmcclure/Documents/GLEE2/output/global_reservoir_emissions.csv"),
                     append = T,
                     row.names = F,
                     col.names = !file.exists("/Users/ryanmcclure/Documents/GLEE2/output/global_reservoir_emissions.csv"))


}

s = Sys.time()

no_cores <- detectCores()
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
foreach(i=country) %dopar% calc_global_flux(i)

e <- Sys.time()
t=e-s
print(t)



# 
# 
# # Aggregated Data For a boxplot
# e <- as.data.frame(do.call(rbind, out_jags)) %>%
#   group_by(date, hylak_id, lat, lon, total_km2) %>%
#   summarize(mean1_4 = mean(ebu_JAGS_FO_1_4),
#             mean_50 = mean(ebu_JAGS_FO_50),
#             mean = mean(ebu_JAGS_base),
#             mean_NLMS = mean(Cal_no_JAGS),
#             mean_Aben = mean(Aben_no_uncertainty)) %>%
#   mutate(mean1_4=ifelse(mean1_4<=0,1,mean1_4),
#          mean_50=ifelse(mean_50<=0,1,mean_50),
#          mean=ifelse(mean<=0,1,mean),
#          mean_NLMS=ifelse(mean_NLMS<=0,1,mean_NLMS),
#          mean_Aben=ifelse(mean<=0,1,mean_Aben)) %>%
#   # Convert to Kg CH4 per year
#   mutate(mean1_4=mean1_4*total_km2*1000000*0.4*365*0.000001,
#          mean_50=mean_50*total_km2*1000000*0.4*365*0.000001,
#          mean=mean*total_km2*1000000*0.4*365*0.000001,
#          mean_NLMS=mean_NLMS*total_km2*1000000*0.4*365*0.000001,
#          mean_Aben=mean_Aben*total_km2*1000000*0.4*365*0.000001) %>%
#   ungroup(.) %>%
#   select(date, mean1_4, mean_50, mean, mean_Aben) %>%
#   reshape2::melt(., id = c("date")) %>%
#   mutate(value = log(value)) %>%
#   filter(!is.infinite(value))
# 
# x[(x == 0)] <- 1
# 
# n <- boxplot(value ~ variable, data = e, 
#              outline = F, 
#              names=c('1-4 Samples', '50+ Samples', 'Sample Error in Model', 'Johnson Paper'),
#              ylab = "log(Kg CH4 per year)",
#              xlab = "Model calibration scenario", 
#              main = "Ebullition flux: all waterbody types in Germany")
# 
# ratio <- tibble::tibble(a = c(diff(n$stats[c(2,4), 1])*0.069305871,diff(n$stats[c(2,4), 1])*0.726284013,diff(n$stats[c(2,4), 1])*0.204410116),
#                         b = c(diff(n$stats[c(2,4), 2])*0.064645843,diff(n$stats[c(2,4), 2])*0.009609274,diff(n$stats[c(2,4), 2])*0.925744883),
#                         c = c(diff(n$stats[c(2,4), 3])*0.07128828,diff(n$stats[c(2,4), 3])*0.92871172, NA))
# # the final plot
# 
# rect(xleft = c(0.6) + seq_along(n$n[1])-1, xright = 1.4 + seq_along(n$n[1])-1, ybottom = n$stats[2, 1]+ratio$a[1]+ratio$a[2]+ratio$a[3],ytop = n$stats[2, 1], col = "#56B4E9")
# rect(xleft = c(0.6) + seq_along(n$n[1])-1, xright = 1.4 + seq_along(n$n[1])-1, ybottom = n$stats[2, 1], ytop = n$stats[2, 1]+ratio$a[1], col = "#D55E00")
# rect(xleft = c(0.6) + seq_along(n$n[1])-1, xright = 1.4 + seq_along(n$n[1])-1, ybottom = n$stats[2, 1]+ratio$a[1], ytop = n$stats[2, 1]+ratio$a[2], col = "#CC79A7")
# 
# 
# rect(xleft = c(0.6) + seq_along(n$n[2]), xright = 1.4 + seq_along(n$n[2]), ybottom = n$stats[2, 2]+ratio$b[1]+ratio$b[2]+ratio$b[3],ytop = n$stats[2, 2], col = "#56B4E9")
# rect(xleft = c(0.6) + seq_along(n$n[2]), xright = 1.4 + seq_along(n$n[2]), ybottom = n$stats[2, 2], ytop = n$stats[2, 2]+ratio$b[1], col = "#D55E00")
# rect(xleft = c(0.6) + seq_along(n$n[2]), xright = 1.4 + seq_along(n$n[2]), ybottom = n$stats[2, 2]+ratio$b[1], ytop = n$stats[2, 2]+ratio$b[2]+ratio$b[1], col = "#CC79A7")
# 
# rect(xleft = c(0.6) + seq_along(n$n[3])+1, xright = 1.4 + seq_along(n$n[3])+1, ybottom = n$stats[2, 3]+ratio$c[1]+ratio$c[2]+ratio$c[1],ytop = n$stats[2, 3], col = "#56B4E9")
# rect(xleft = c(0.6) + seq_along(n$n[3])+1, xright = 1.4 + seq_along(n$n[3])+1, ybottom = n$stats[2, 3], ytop = n$stats[2, 3]+ratio$c[1], col = "#D55E00")
# 
