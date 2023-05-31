# =======================================================================
#------------------------------------------------------------------------

#### Libraries #### 
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(vroom, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)

Aben_model_paper <- function(temp){
  est = 100 * 1.1 ^ (temp-20)
  return(est)
}

Aben_model_NMLS <- function(temp){
  est = (rnorm(1,77.07198, sd = 8.55152) * (rnorm(1,1.04613, sd = 0.01832) ^ (temp-20))) + rnorm(1,0, sd = 120)
  return(est)
}

ebu_arhennius <- function(temp){
  est = rnorm(1,98.44878, sd = 16.16218) * rnorm(1,1.145583, sd = 0.02777339) ^ (temp-20) + rnorm(1,0, sd = 198.1507)
  return(est)
}

ebu_arhennius_FO_1_4 <- function(temp){
  est = rnorm(1,95.89238, sd = 24.677516) * rnorm(1,1.175983, sd = 0.05976376) ^ (temp-20) + rnorm(1,0, sd = 224.1444660) + rnorm(1,0, sd = 56.98130)
  return(est)
}

ebu_arhennius_FO_50 <- function(temp){
  est = rnorm(1,97.12485, sd = 8.341771) * rnorm(1,1.171232, sd = 0.01409348) ^ (temp-20) + rnorm(1,0, sd = 0.9928644) + rnorm(1,0, sd = 99.93237)
  return(est)
}

file1 <- "/Users/ryanmcclure/Documents/GLEE-v1.0/data/GLEE-upscale-countries/Brazil.csv"

file2 <- "/Users/ryanmcclure/Documents/GLEE-v1.0/data/GLEE-upscale-countries/Germany.csv"

output <- "/central/groups/carnegie_poc/rmcclure/project-GLEE/global_ebullition_emission.csv"

d <- vroom::vroom(paste0(file2), col_names = F) %>%
  filter(X1 >= 2002) %>%
  mutate(X9 = ifelse(is.na(X9),0,X9)) %>%
  na.omit(.)
  
#calculate_emission <- function(country){

out_jags <- list()

iteration <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21)

for(i in 1:length(iteration)){
  
  d1 <- d %>%
    rename(year = X1, month= X2,
           hylak_id = X3, lat = X4,
           lon = X5, country = X6,
           total_km2 = X7, mean_temp_k = X8,
           ice_cover_mean = X9) %>%
    group_by(hylak_id, lat, lon, country, total_km2) %>%
    mutate(ebu_JAGS_FO_1_4 = ifelse(ice_cover_mean > 50, 0, mapply(ebu_arhennius_FO_1_4, mean_temp_k-273.15)),
           ebu_JAGS_FO_50 = ifelse(ice_cover_mean > 50, 0, mapply(ebu_arhennius_FO_50, mean_temp_k-273.15)),
           ebu_JAGS_base = ifelse(ice_cover_mean > 50, 0, mapply(ebu_arhennius, mean_temp_k-273.15)),
           Cal_no_JAGS = ifelse(ice_cover_mean > 50, 0, mapply(Aben_model_NMLS, mean_temp_k-273.15)),
           Aben_no_uncertainty = ifelse(ice_cover_mean > 50, 0, mapply(Aben_model_paper, mean_temp_k-273.15)),
           'date' = lubridate::make_date(year = year, month = month),
           iteration = iteration[i])
  
  out_jags[[i]] <- d1
  
}

# Raw Data
# b <- as.data.frame(do.call(rbind, out_jags)) %>%
#   mutate(ebu_JAGS_FO_50=ifelse(ebu_JAGS_FO_50<=0,0,ebu_JAGS_FO_50),
#          ebu_JAGS_FO_1_4=ifelse(ebu_JAGS_FO_1_4<=0,0,ebu_JAGS_FO_1_4),
#          ebu_JAGS_base=ifelse(ebu_JAGS_base<=0,0,ebu_JAGS_base),
#          Cal_no_JAGS=ifelse(mean_NLMS<=0,0,Cal_no_JAGS),
#          Aben_no_uncertainty=ifelse(mean<=0,0,Aben_no_uncertainty))


# Aggregated Data
e <- as.data.frame(do.call(rbind, out_jags)) %>%
  group_by(date, hylak_id, lat, lon, total_km2) %>%
  summarize(mean1_4 = mean(ebu_JAGS_FO_1_4),
            mean_50 = mean(ebu_JAGS_FO_50),
            mean = mean(ebu_JAGS_base),
            mean_NLMS = mean(Cal_no_JAGS),
            mean_Aben = mean(Aben_no_uncertainty)) %>%
  mutate(mean1_4=ifelse(mean1_4<=0,1,mean1_4),
         mean_50=ifelse(mean_50<=0,1,mean_50),
         mean=ifelse(mean<=0,1,mean),
         mean_NLMS=ifelse(mean_NLMS<=0,1,mean_NLMS),
         mean_Aben=ifelse(mean<=0,1,mean_Aben)) %>%
  # Convert to Kg CH4 per year
  mutate(mean1_4=mean1_4*total_km2*1000000*0.4*365*0.000001,
         mean_50=mean_50*total_km2*1000000*0.4*365*0.000001,
         mean=mean*total_km2*1000000*0.4*365*0.000001,
         mean_NLMS=mean_NLMS*total_km2*1000000*0.4*365*0.000001,
         mean_Aben=mean_Aben*total_km2*1000000*0.4*365*0.000001) %>%
  ungroup(.) %>%
  select(date, mean1_4, mean_50, mean, mean_Aben) %>%
  reshape2::melt(., id = c("date")) %>%
  mutate(value = log(value)) %>%
  filter(!is.infinite(value))

x[(x == 0)] <- 1

n <- boxplot(value ~ variable, data = e, 
             outline = F, 
             names=c('1-4 Samples', '50+ Samples', 'Sample Error in Model', 'Johnson Paper'),
             ylab = "log(Kg CH4 per year)",
             xlab = "Model calibration scenario", 
             main = "Ebullition flux: all waterbody types in Germany")

ratio <- tibble::tibble(a = c(diff(n$stats[c(2,4), 1])*0.069305871,diff(n$stats[c(2,4), 1])*0.726284013,diff(n$stats[c(2,4), 1])*0.204410116),
                        b = c(diff(n$stats[c(2,4), 2])*0.064645843,diff(n$stats[c(2,4), 2])*0.009609274,diff(n$stats[c(2,4), 2])*0.925744883),
                        c = c(diff(n$stats[c(2,4), 3])*0.07128828,diff(n$stats[c(2,4), 3])*0.92871172, NA))
# the final plot

rect(xleft = c(0.6) + seq_along(n$n[1])-1, xright = 1.4 + seq_along(n$n[1])-1, ybottom = n$stats[2, 1]+ratio$a[1]+ratio$a[2]+ratio$a[3],ytop = n$stats[2, 1], col = "#56B4E9")
rect(xleft = c(0.6) + seq_along(n$n[1])-1, xright = 1.4 + seq_along(n$n[1])-1, ybottom = n$stats[2, 1], ytop = n$stats[2, 1]+ratio$a[1], col = "#D55E00")
rect(xleft = c(0.6) + seq_along(n$n[1])-1, xright = 1.4 + seq_along(n$n[1])-1, ybottom = n$stats[2, 1]+ratio$a[1], ytop = n$stats[2, 1]+ratio$a[2], col = "#CC79A7")


rect(xleft = c(0.6) + seq_along(n$n[2]), xright = 1.4 + seq_along(n$n[2]), ybottom = n$stats[2, 2]+ratio$b[1]+ratio$b[2]+ratio$b[3],ytop = n$stats[2, 2], col = "#56B4E9")
rect(xleft = c(0.6) + seq_along(n$n[2]), xright = 1.4 + seq_along(n$n[2]), ybottom = n$stats[2, 2], ytop = n$stats[2, 2]+ratio$b[1], col = "#D55E00")
rect(xleft = c(0.6) + seq_along(n$n[2]), xright = 1.4 + seq_along(n$n[2]), ybottom = n$stats[2, 2]+ratio$b[1], ytop = n$stats[2, 2]+ratio$b[2]+ratio$b[1], col = "#CC79A7")

rect(xleft = c(0.6) + seq_along(n$n[3])+1, xright = 1.4 + seq_along(n$n[3])+1, ybottom = n$stats[2, 3]+ratio$c[1]+ratio$c[2]+ratio$c[1],ytop = n$stats[2, 3], col = "#56B4E9")
rect(xleft = c(0.6) + seq_along(n$n[3])+1, xright = 1.4 + seq_along(n$n[3])+1, ybottom = n$stats[2, 3], ytop = n$stats[2, 3]+ratio$c[1], col = "#D55E00")







mutate(sd_low_1_4 = ifelse(mean_1_4-sd_1_4 < 0,0,sd_1_4),
         sd_low_1_4 = ifelse(sd_low_1_4 == 0, mean_1_4, sd_low_1_4),
         sd_low_50 = ifelse(mean_50-sd_50 < 0,0,sd_50),
         sd_low_50 = ifelse(sd_low_50 == 0, mean_50, sd_low_50),
         sd_low = ifelse(mean-sd < 0,0,sd),
         sd_low = ifelse(sd_low == 0, mean, sd_low)) %>%
  reshape2::melt(., id = c("date", "hylak_id", "lat", "lon", "country")) %>%
    utils::write.table(., file = paste0("/Users/ryanmcclure/Documents/GLEE2/output/global_emissions.csv"),
                       append = T,
                       row.names = F,
                       col.names = !file.exists("/Users/ryanmcclure/Documents/GLEE2/output/global_emissions.csv"))
  


s = Sys.time()

no_cores <- detectCores()
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
foreach(i=country) %dopar% calculate_emission(i)

e <- Sys.time()
t=e-s
print(t)



# p1 <- JAGS_output_summary %>% filter(variable == "mean_1_4") %>%
# ggplot(., aes(x = date, y = value))+
#   geom_line(lwd = 0.5, alpha = 1)+
#   geom_ribbon(aes(x = date, y = value, ymin = value-sd_low_1_4, ymax = value+sd_1_4), alpha = 0.2, lwd = 0)+
#   geom_point(data = ebu_obs_compare, aes(date, ch4_ebu), pch = 21, color = "grey", 
#              fill = "red", inherit.aes = F)+
#   coord_cartesian(ylim = c(-10,400))+
#   labs(title = "Assuming Low Observation Certainty (i.e. 1-4 total samples observed)")+
#   theme_classic()

# p1 <- JAGS_output_summary %>% filter(variable == "mean_1_4_LAC") %>%
#   mutate(value = ifelse(value < 0, 0, value * 0.000001),
#          sd_1_4_LAC = sd_1_4_LAC * 0.000001,
#          sd_low_1_4_LAC = sd_low_1_4_LAC* 0.000001) %>%
#   ggplot(., aes(x = date, y = value))+
#   geom_line(lwd = 0.5, alpha = 1)+
#   geom_ribbon(aes(x = date, y = value, ymin = value-sd_low_1_4_LAC, ymax = value+sd_1_4_LAC), alpha = 0.2, lwd = 0)+
#   labs(title = "Includes Low Observation Certainty (1-4) and Model Process Error")+
#   ylab("kg CH4/day")+
#   ylim(0,350)+
#   theme_classic()

# p2 <- JAGS_output_summary %>% filter(variable == "mean_50") %>%
#   ggplot(., aes(x = date, y = value))+
#   geom_line(lwd = 0.5, alpha = 1)+
#   geom_ribbon(aes(x = date, y = value, ymin = value-sd_low_50, ymax = value+sd_50), alpha = 0.2, lwd = 0)+
#   geom_point(data = ebu_obs_compare, aes(date, ch4_ebu), pch = 21, color = "grey", 
#              fill = "red", inherit.aes = F)+
#   labs(title = "Assuming High Observation Certainty (i.e. 50+ samples) and Includes Model Process Error")+
#   ylim(0,350)+
#   theme_classic()

# p2 <- JAGS_output_summary %>% filter(variable == "mean_50_LAC") %>%
#   mutate(value = ifelse(value < 0, 0, value * 0.000001),
#          sd_50_LAC = sd_50_LAC * 0.000001,
#          sd_low_50_LAC = sd_low_50_LAC* 0.000001) %>%
#   ggplot(., aes(x = date, y = value))+
#   geom_line(lwd = 0.5, alpha = 1)+
#   geom_ribbon(aes(x = date, y = value, ymin = value-sd_low_50_LAC, ymax = value+sd_50_LAC), alpha = 0.2, lwd = 0)+
#   labs(title = "Includes High Observation Certainty (50+) and Model Process Error")+
#   ylab("kg CH4/day")+
#   ylim(0,350)+
#   theme_classic()











# JAGS_output_summary %>% filter(variable == "mean_50_LAC") %>%
#   summarise(mean_value = sum(value * 0.000001*30)/13.74,
#          sd_50_LAC_high = sum(sd_50_LAC * 0.000001*30)/13.74,
#          sd_low_50_LAC = sum((value-sd_low_50_LAC)* 0.000001*30)/13.74)
# 
# JAGS_output_summary %>% filter(variable == "mean_1_4_LAC") %>%
#   summarise(mean_value = sum(value * 0.000001*30)/13.74,
#             sd_1_4_LAC_high = sum(sd_1_4_LAC * 0.000001*30)/13.74,
#             sd_low_1_4_LAC = sum((value-sd_low_1_4_LAC)* 0.000001*30)/13.74)



