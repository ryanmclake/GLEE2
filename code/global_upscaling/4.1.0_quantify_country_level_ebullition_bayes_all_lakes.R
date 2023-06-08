# Clean your environment
rm(list=ls())

#Clean your memory
gc()

# Set seed for repeat measurements
set.seed(3818)

# =======================================================================
#------------------------------------------------------------------------

#### Libraries #### 
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(vroom, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)

Johnson_model_ebullition <- function(temp){
  est = 100 * 1.1 ^ (temp-20)
  return(est)
}

Johnson_model_diffusion <- function(temp){
  est = 0.023*exp(0.124*temp)
  return(est)
}

arrhenius_ebu_posterior_prediction <- function(E20, omega, temp, Q){
  est = (E20 * omega ^ (temp-20)) + rnorm(1000, 0, sd = Q)
  return(est)
}

arrhenius_diff_posterior_prediction <- function(A, a, temp, Q){
  est = (A*exp(a*(temp))) + rnorm(1000, 0, sd = Q)
  return(est)
}

ebu_coefficients <- read_csv("./output/raw_ebullition_coefficients.csv")
diff_coefficients <- read_csv("./output/raw_diffusion_coefficients.csv")

file <- "/Users/ryanmcclure/Documents/GLEE2.1/glcp_countries/United States of America.csv"

d <- vroom::vroom(paste0(file), col_names = F) %>%
  mutate(temp_for_model_C = X10-273.15)
  
temp_data <- c(d$temp_for_model_C)

parms <- ebu_coefficients %>%
  filter(waterbody_type == "lake") %>%
  sample_n(., 100, replace=TRUE)

predict_output <- list()

for(s in 1:length(temp_data)){
  
  prediction <- arrhenius_ebu_posterior_prediction(temp = temp_data[s],
                                                   E20 = parms$E20,
                                                   omega = parms$omega,
                                                   Q = parms$sd.pro)
  predict_output[[s]] <- prediction
}

prediction_output = as.data.frame(do.call(rbind, predict_output))

prediction_output_long <- prediction_output %>% t(.) %>% reshape2::melt(.) %>%
  group_by(Var2) %>%
  summarize(mean = mean(value),
            sd = sd(value),
            var = var(value))

#calculate_emission <- function(country){

out_jags <- list()

iteration <- seq(1:50)

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

# Raw Data to plot with shadows
e2 <- as.data.frame(do.call(rbind, out_jags)) %>%
  group_by(date) %>%
  summarize(mean_1_4 = mean(ebu_JAGS_FO_1_4),
            mean_50 = mean(ebu_JAGS_FO_50),
            mean = mean(ebu_JAGS_base),
            mean_Aben = mean(Aben_no_uncertainty),
            sd_1_4 = sd(ebu_JAGS_FO_1_4),
            sd_50 = sd(ebu_JAGS_FO_50),
            sd = sd(ebu_JAGS_base),
            sd_Aben = sd(Aben_no_uncertainty)) %>%
  mutate(mean_1_4=ifelse(mean_1_4<=0,1,mean_1_4),
         mean_50=ifelse(mean_50<=0,1,mean_50),
         mean=ifelse(mean<=0,1,mean),
         mean_Aben=ifelse(mean<=0,1,mean_Aben)) %>%
  mutate(sd_low_1_4 = ifelse(mean_1_4-sd_1_4 < 0,0,sd_1_4),
         sd_low_1_4 = ifelse(sd_low_1_4 == 0, mean_1_4, sd_low_1_4),
         sd_low_50 = ifelse(mean_50-sd_50 < 0,0,sd_50),
         sd_low_50 = ifelse(sd_low_50 == 0, mean_50, sd_low_50),
         sd_low = ifelse(mean-sd < 0,0,sd),
         sd_low = ifelse(sd_low == 0, mean, sd_low),
         sd_Aben_low = ifelse(mean_Aben-sd_Aben < 0,0,sd_Aben),
         sd_Aben_low = ifelse(sd_low == 0, mean_Aben, sd_Aben_low)) %>%
        reshape2::melt(., id = c("date","sd_low_1_4","sd_low_50","sd_low","sd_Aben_low",
                                 "sd_1_4","sd_50","sd","sd_Aben"))

 p1 <- e2 %>% filter(variable == "mean_1_4") %>%
   select(date, value, sd_low_1_4, sd_1_4) %>%
   ggplot(., aes(x = date, y = value))+
   geom_line(lwd = 0.5, alpha = 1)+
   geom_ribbon(aes(x = date, y = value, ymin = value-sd_low_1_4, ymax = value+sd_1_4), alpha = 0.2, lwd = 0)+
   #geom_point(data = ebu_obs_compare, aes(date, ch4_ebu), pch = 21, color = "grey",
   #            fill = "red", inherit.aes = F)+
   coord_cartesian(ylim = c(-10,400))+
   labs(title = "Assuming High Measurement Unertainty (i.e. 1-4 total samples observed)")+
   theme_classic()

 p2 <- e2 %>% filter(variable == "mean") %>%
   select(date, value, sd_low, sd) %>%
   ggplot(., aes(x = date, y = value))+
   geom_line(lwd = 0.5, alpha = 1)+
   geom_ribbon(aes(x = date, y = value, ymin = value-sd_low, ymax = value+sd), alpha = 0.2, lwd = 0)+
   #geom_point(data = ebu_obs_compare, aes(date, ch4_ebu), pch = 21, color = "grey",
   #            fill = "red", inherit.aes = F)+
   coord_cartesian(ylim = c(-10,400))+
   labs(title = "Intregrate Measurement Error into overall Model Uncertainty")+
   theme_classic()
 
 p3 <- e2 %>% filter(variable == "mean_50") %>%
   select(date, value, sd_low_50, sd_50) %>%
   ggplot(., aes(x = date, y = value))+
   geom_line(lwd = 0.5, alpha = 1)+
   geom_ribbon(aes(x = date, y = value, ymin = value-sd_low_50, ymax = value+sd_50), alpha = 0.2, lwd = 0)+
   #geom_point(data = ebu_obs_compare, aes(date, ch4_ebu), pch = 21, color = "grey",
   #            fill = "red", inherit.aes = F)+
   coord_cartesian(ylim = c(-10,400))+
   labs(title = "Assuming Low Measurement Unertainty (i.e. 50+ total samples observed)")+
   theme_classic()

 p4 <- e2 %>% filter(variable == "mean_Aben") %>%
   select(date, value, sd_Aben_low, sd_Aben) %>%
   ggplot(., aes(x = date, y = value))+
   geom_line(lwd = 0.5, alpha = 1)+
   geom_ribbon(aes(x = date, y = value, ymin = value-sd_Aben_low, ymax = value+sd_Aben), alpha = 0.2, lwd = 0)+
   #geom_point(data = ebu_obs_compare, aes(date, ch4_ebu), pch = 21, color = "grey",
   #            fill = "red", inherit.aes = F)+
   coord_cartesian(ylim = c(-10,400))+
   labs(title = "Using Aben Model")+
   theme_classic()
 
 library(patchwork)
 
 plot <- (p2+p4)/(p1+p3)

# Aggregated Data For a boxplot
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

