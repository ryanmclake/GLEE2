library(minpack.lm)
library(hydroGOF)

set.seed(1098)

#How does observation, coefficient, and model variability contribute to total global emissions uncertainty?
### EXPLORE ERROR IN THE OBSERVATIONS FOR THE Q10 TEMP MODELS ###
# Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)

ebu_base_temp <- error_ebu %>%
  ungroup(.) %>%
  select(ch4_ebu, sd_time, waterbody_type, sd_space, temp_for_model_C, T_OC, mean_sw_wm2) %>%
  mutate(ebu_time_low = ifelse(ch4_ebu-sd_time < 0, 0.01, ch4_ebu-sd_time)) %>%
  mutate(ebu_time_high = ifelse(ch4_ebu+sd_time < 0, 0.01, ch4_ebu+sd_time)) %>%
  mutate(ebu_space_low = ifelse(ch4_ebu-sd_space < 0, 0.01, ch4_ebu-sd_space)) %>%
  mutate(ebu_space_high = ifelse(ch4_ebu+sd_space < 0, 0.01, ch4_ebu+sd_space)) %>%
  na.omit(.)%>%
  filter(waterbody_type == "lake")

mean_ebu_MFOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
                      start = list(A = 100, a = 1.1),
                      data = ebu_base_temp,
                      control = nls.lm.control(maxiter=1000))


summary(mean_ebu_MFOA) # get MSE value

eval_models <- as.data.frame(cbind(ebu_base_temp$ch4_ebu, predict(mean_ebu_MFOA)))

eval_FOA_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_main_a <- rbind(mean_ebu_MFOA$m$getPars())


coefficients_low_a <- cbind(mean_ebu_MFOA$m$getPars()[1]-summary(mean_ebu_MFOA)$coefficients[,2][1],
                          mean_ebu_MFOA$m$getPars()[2]-summary(mean_ebu_MFOA)$coefficients[,2][2])

coefficients_high_a <- cbind(mean_ebu_MFOA$m$getPars()[1]+summary(mean_ebu_MFOA)$coefficients[,2][1],
                           mean_ebu_MFOA$m$getPars()[2]+summary(mean_ebu_MFOA)$coefficients[,2][2])

FOA_mean <- function(x){
  predicted_ebu_rate = (coefficients_main_a[1,1] * coefficients_main_a[1,2]^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_low <- function(x){
  predicted_ebu_rate = (coefficients_low_a[1,1] * coefficients_low_a[1,2]^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_high <- function(x){
  predicted_ebu_rate = (coefficients_high_a[1,1] * coefficients_high_a[1,2]^(x-20))
  return(predicted_ebu_rate)
}


mean_ebu_SOA_temp_soil_threshold = nlsLM(ch4_ebu ~ A * exp((a * temp_for_model_C)*(T_OC/k+T_OC)),
                                         start = list(A = 100, a = 0.035, k = 40),
                                         data = ebu_base_temp,
                                         control = nls.lm.control(maxiter=1024))


summary(mean_ebu_SOA_temp_soil_threshold) # get MSE value

eval_models <- as.data.frame(cbind(ebu_base_temp$ch4_ebu, predict(mean_ebu_SOA_temp_soil_threshold)))

eval_SOA_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_main_b <- rbind(mean_ebu_SOA_temp_soil_threshold$m$getPars())


coefficients_low_b <- cbind(mean_ebu_SOA_temp_soil_threshold$m$getPars()[1]-summary(mean_ebu_SOA_temp_soil_threshold)$coefficients[,2][1],
                          mean_ebu_SOA_temp_soil_threshold$m$getPars()[2]-summary(mean_ebu_SOA_temp_soil_threshold)$coefficients[,2][2],
                          mean_ebu_SOA_temp_soil_threshold$m$getPars()[3]-summary(mean_ebu_SOA_temp_soil_threshold)$coefficients[,2][3])

coefficients_high_b <- cbind(mean_ebu_SOA_temp_soil_threshold$m$getPars()[1]+summary(mean_ebu_SOA_temp_soil_threshold)$coefficients[,2][1],
                           mean_ebu_SOA_temp_soil_threshold$m$getPars()[2]+summary(mean_ebu_SOA_temp_soil_threshold)$coefficients[,2][2],
                           mean_ebu_SOA_temp_soil_threshold$m$getPars()[3]+summary(mean_ebu_SOA_temp_soil_threshold)$coefficients[,2][3])

SOA_coefficient_mean <- function(x, p){
  predicted_ebu_rate = coefficients_main_b[1,1] * exp((coefficients_main_b[1,2]*x) * (p/coefficients_main_b[1,3]+p))
  return(predicted_ebu_rate)
}

SOA_coefficient_low <- function(x, p){
  predicted_ebu_rate = coefficients_low_b[1,1] * exp((coefficients_low_b[1,2]*x) * (p/coefficients_low_b[1,3]+p))
  return(predicted_ebu_rate)
}

SOA_coefficient_high <- function(x, p){
  predicted_ebu_rate = coefficients_high_b[1,1] * exp((coefficients_high_b[1,2]*x) * (p/coefficients_high_b[1,3]+p))
  return(predicted_ebu_rate)
}



LM_mean_ebu = glm(ch4_ebu ~ mean_sw_wm2, data = ebu_base_temp)

summary(LM_mean_ebu) # get MSE value

eval_models <- as.data.frame(cbind(ebu_base_temp$ch4_ebu, predict(LM_mean_ebu)))

eval_LM_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_mean_c <- rbind(LM_mean_ebu$coefficients)

coefficients_low_c <- cbind(LM_mean_ebu$coefficients[1]-summary(LM_mean_ebu)$coefficients[,2][1],
                          LM_mean_ebu$coefficients[2]-summary(LM_mean_ebu)$coefficients[,2][2])

coefficients_high_c <- cbind(LM_mean_ebu$coefficients[1]+summary(LM_mean_ebu)$coefficients[,2][1],
                          LM_mean_ebu$coefficients[2]+summary(LM_mean_ebu)$coefficients[,2][2])


LM_coefficient_mean <- function(x, p){
  predicted_ebu_rate = coefficients_mean_c[1,1] + p*coefficients_mean_c[1,2]
  return(predicted_ebu_rate)
}

LM_coefficient_low <- function(x, p){
  predicted_ebu_rate = coefficients_low_c[1,1] + p*coefficients_low_c[1,2]
  return(predicted_ebu_rate)
}

LM_coefficient_high <- function(x, p){
  predicted_ebu_rate = coefficients_high_c[1,1] + p*coefficients_high_c[1,2]
  return(predicted_ebu_rate)
}






ebu_base_temp <- error_ebu %>%
  ungroup(.) %>%
  select(ch4_ebu, sd_time, waterbody_type, sd_space, temp_for_model_C, T_OC, mean_sw_wm2) %>%
  mutate(ebu_time_low = ifelse(ch4_ebu-sd_time < 0, 0.01, ch4_ebu-sd_time)) %>%
  mutate(ebu_time_high = ifelse(ch4_ebu+sd_time < 0, 0.01, ch4_ebu+sd_time)) %>%
  mutate(ebu_space_low = ifelse(ch4_ebu-sd_space < 0, 0.01, ch4_ebu-sd_space)) %>%
  mutate(ebu_space_high = ifelse(ch4_ebu+sd_space < 0, 0.01, ch4_ebu+sd_space)) %>%
  na.omit(.)%>%
  filter(waterbody_type == "reservoir")

mean_ebu_MFOA_res = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
                      start = list(A = 100, a = 1.1),
                      data = ebu_base_temp,
                      control = nls.lm.control(maxiter=1000))


summary(mean_ebu_MFOA_res) # get MSE value

eval_models <- as.data.frame(cbind(ebu_base_temp$ch4_ebu, predict(mean_ebu_MFOA_res)))

eval_FOA_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_main_d <- rbind(mean_ebu_MFOA_res$m$getPars())


coefficients_low_d <- cbind(mean_ebu_MFOA_res$m$getPars()[1]-summary(mean_ebu_MFOA_res)$coefficients[,2][1],
                          mean_ebu_MFOA_res$m$getPars()[2]-summary(mean_ebu_MFOA_res)$coefficients[,2][2])

coefficients_high_d <- cbind(mean_ebu_MFOA_res$m$getPars()[1]+summary(mean_ebu_MFOA_res)$coefficients[,2][1],
                           mean_ebu_MFOA_res$m$getPars()[2]+summary(mean_ebu_MFOA_res)$coefficients[,2][2])

FOA_coefficient_mean_res <- function(x){
  predicted_ebu_rate = (coefficients_main_d[1,1] * coefficients_main_d[1,2]^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_low_res <- function(x){
  predicted_ebu_rate = (coefficients_low_d[1,1] * coefficients_low_d[1,2]^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_high_res <- function(x){
  predicted_ebu_rate = (coefficients_high_d[1,1] * coefficients_high_d[1,2]^(x-20))
  return(predicted_ebu_rate)
}


mean_ebu_SOA_temp_soil_threshold_res = nlsLM(ch4_ebu ~ A * exp((a * temp_for_model_C)*(T_OC/k+T_OC)),
                                         start = list(A = 100, a = 0.035, k = 40),
                                         data = ebu_base_temp,
                                         control = nls.lm.control(maxiter=1024))


summary(mean_ebu_SOA_temp_soil_threshold_res) # get MSE value

eval_models <- as.data.frame(cbind(ebu_base_temp$ch4_ebu, predict(mean_ebu_SOA_temp_soil_threshold_res)))

eval_SOA_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_main_e <- rbind(mean_ebu_SOA_temp_soil_threshold_res$m$getPars())


coefficients_low_e <- cbind(mean_ebu_SOA_temp_soil_threshold_res$m$getPars()[1]-summary(mean_ebu_SOA_temp_soil_threshold_res)$coefficients[,2][1],
                          mean_ebu_SOA_temp_soil_threshold_res$m$getPars()[2]-summary(mean_ebu_SOA_temp_soil_threshold_res)$coefficients[,2][2],
                          mean_ebu_SOA_temp_soil_threshold_res$m$getPars()[3]-summary(mean_ebu_SOA_temp_soil_threshold_res)$coefficients[,2][3])

coefficients_high_e <- cbind(mean_ebu_SOA_temp_soil_threshold_res$m$getPars()[1]+summary(mean_ebu_SOA_temp_soil_threshold_res)$coefficients[,2][1],
                           mean_ebu_SOA_temp_soil_threshold_res$m$getPars()[2]+summary(mean_ebu_SOA_temp_soil_threshold_res)$coefficients[,2][2],
                           mean_ebu_SOA_temp_soil_threshold_res$m$getPars()[3]+summary(mean_ebu_SOA_temp_soil_threshold_res)$coefficients[,2][3])

SOA_coefficient_mean_res <- function(x, p){
  predicted_ebu_rate = coefficients_main_e[1,1] * exp((coefficients_main_e[1,2]*x) * (p/coefficients_main_e[1,3]+p))
  return(predicted_ebu_rate)
}

SOA_coefficient_low_res <- function(x, p){
  predicted_ebu_rate = coefficients_low_e[1,1] * exp((coefficients_low_e[1,2]*x) * (p/coefficients_low_e[1,3]+p))
  return(predicted_ebu_rate)
}

SOA_coefficient_high_res <- function(x, p){
  predicted_ebu_rate = coefficients_high_e[1,1] * exp((coefficients_high_e[1,2]*x) * (p/coefficients_high_e[1,3]+p))
  return(predicted_ebu_rate)
}



LM_mean_ebu_res = glm(ch4_ebu ~ mean_sw_wm2, data = ebu_base_temp)

summary(LM_mean_ebu_res) # get MSE value

eval_models <- as.data.frame(cbind(ebu_base_temp$ch4_ebu, predict(LM_mean_ebu_res)))

eval_LM_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_mean_f <- rbind(LM_mean_ebu_res$coefficients)

coefficients_low_f <- cbind(LM_mean_ebu_res$coefficients[1]-summary(LM_mean_ebu_res)$coefficients[,2][1],
                          LM_mean_ebu_res$coefficients[2]-summary(LM_mean_ebu_res)$coefficients[,2][2])

coefficients_high_f <- cbind(LM_mean_ebu_res$coefficients[1]+summary(LM_mean_ebu_res)$coefficients[,2][1],
                           LM_mean_ebu_res$coefficients[2]+summary(LM_mean_ebu_res)$coefficients[,2][2])


LM_coefficient_mean_res <- function(x, p){
  predicted_ebu_rate = coefficients_mean_f[1,1] + p*coefficients_mean_f[1,2]
  return(predicted_ebu_rate)
}

LM_coefficient_low_res <- function(x, p){
  predicted_ebu_rate = coefficients_low_f[1,1] + p*coefficients_low_f[1,2]
  return(predicted_ebu_rate)
}

LM_coefficient_high_res <- function(x, p){
  predicted_ebu_rate = coefficients_high_f[1,1] + p*coefficients_high_f[1,2]
  return(predicted_ebu_rate)
}





country <- c(list.files("/Users/ryanmcclure/Documents/GLEE2.1/data/countries"))

calc_global_variance <- function(country){
  
  d <- vroom::vroom(paste0("/Users/ryanmcclure/Documents/GLEE2.1/data/United States of America.csv"), col_names = F, delim = ",") %>%
    rename(year = X1, month= X2,
           hylak_id = X3, lat = X4,
           lon = X5, continent = X6, country = X7,
           total_km2 = X8, waterbody_type = X9,
           mean_temp_k = X10, ice_cover_mean = X11) %>%
    mutate(ice_cover_mean = ifelse(is.na(ice_cover_mean),0,ice_cover_mean)) %>%
    na.omit(.)
  
  d_model <- d %>%
    group_by(hylak_id, lat, lon, total_km2) %>%
    mutate(FOA_mean_ebu = ifelse(waterbody_type == 1, mapply(FOA_mean, x=mean_temp_k-273.15), mapply(FOA_coefficient_mean_res, x=mean_temp_k-273.15)),
           FOA_coefficient_high_ebu = ifelse(waterbody_type == 1, mapply(FOA_coefficient_high, x=mean_temp_k-273.15), mapply(FOA_coefficient_high_res, x=mean_temp_k-273.15)),
           FOA_coefficient_low_ebu = ifelse(waterbody_type == 1, mapply(FOA_coefficient_low, x=mean_temp_k-273.15), mapply(FOA_coefficient_low_res, x=mean_temp_k-273.15)),
           'date' = lubridate::make_date(year = year, month = month)) %>%
    ungroup(.) %>%
    select(-year, -month) %>%
    utils::write.table(., file = paste0("/Users/ryanmcclure/Documents/GLEE2.1/output/global_ebullition_coefficient_variance.csv"),
                       append = T,
                       row.names = F,
                       col.names = !file.exists("/Users/ryanmcclure/Documents/GLEE2.1/output/global_ebullition_coefficient_variance.csv"))
}


s = Sys.time()

no_cores <- detectCores() - 4
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
foreach(i=country) %dopar% calc_global_variance(i)

e <- Sys.time()
t=e-s
print(t)


measurement <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/global_ebullition_coefficient_variance.csv") %>%
  group_by(lat, lon, hylak_id) %>%
  summarize(FOA_coefficient_high_ebu = median(FOA_coefficient_high_ebu, na.rm = T),
            FOA_coefficient_low_ebu = median(FOA_coefficient_low_ebu, na.rm = T)) %>%
  filter(lon < 60) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")

boxplot_coefficient <- vroom::vroom("/Users/ryanmcclure/Documents/GLEE2.1/output/global_ebullition_coefficient_variance.csv") %>%
  group_by(lat, lon, hylak_id) %>%
  summarize(FOA_coefficient_high_ebu = median(FOA_coefficient_high_ebu, na.rm = T),
            FOA_coefficient_low_ebu = median(FOA_coefficient_low_ebu, na.rm = T)) %>%
  ungroup(.) %>%
  select(hylak_id, FOA_coefficient_high_ebu, FOA_coefficient_low_ebu) %>%
  reshape2::melt(., id = "hylak_id")

ggplot(boxplot_coefficient, aes(x=variable, y=value)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16,
               outlier.size=0.5)

boxes <- rbind(boxplot, boxplot_model, boxplot_coefficient)

ggplot(boxes, aes(x=variable, y=value)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=16,
               outlier.size=0.2)+
  coord_cartesian(ylim = c(0,1200))
