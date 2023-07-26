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

coefficients_main <- rbind(mean_ebu_MFOA$m$getPars())


coefficients_low <- cbind(mean_ebu_MFOA$m$getPars()[1]-summary(mean_ebu_MFOA)$coefficients[,2][1],
                          mean_ebu_MFOA$m$getPars()[2]-summary(mean_ebu_MFOA)$coefficients[,2][2])

coefficients_high <- cbind(mean_ebu_MFOA$m$getPars()[1]+summary(mean_ebu_MFOA)$coefficients[,2][1],
                           mean_ebu_MFOA$m$getPars()[2]+summary(mean_ebu_MFOA)$coefficients[,2][2])

FOA_mean <- function(x){
  predicted_ebu_rate = (coefficients_main[1,1] * coefficients_main[1,2]^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_low <- function(x){
  predicted_ebu_rate = (coefficients_low[1,1] * coefficients_low[1,2]^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_high <- function(x){
  predicted_ebu_rate = (coefficients_high[1,1] * coefficients_high[1,2]^(x-20))
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

coefficients_main <- rbind(mean_ebu_SOA_temp_soil_threshold$m$getPars())


coefficients_low <- cbind(mean_ebu_SOA_temp_soil_threshold$m$getPars()[1]-summary(mean_ebu_SOA_temp_soil_threshold)$coefficients[,2][1],
                          mean_ebu_SOA_temp_soil_threshold$m$getPars()[2]-summary(mean_ebu_SOA_temp_soil_threshold)$coefficients[,2][2],
                          mean_ebu_SOA_temp_soil_threshold$m$getPars()[3]-summary(mean_ebu_SOA_temp_soil_threshold)$coefficients[,2][3])

coefficients_high <- cbind(mean_ebu_SOA_temp_soil_threshold$m$getPars()[1]+summary(mean_ebu_SOA_temp_soil_threshold)$coefficients[,2][1],
                           mean_ebu_SOA_temp_soil_threshold$m$getPars()[2]+summary(mean_ebu_SOA_temp_soil_threshold)$coefficients[,2][2],
                           mean_ebu_SOA_temp_soil_threshold$m$getPars()[3]+summary(mean_ebu_SOA_temp_soil_threshold)$coefficients[,2][3])

SOA_coefficient_mean <- function(x, p){
  predicted_ebu_rate = coefficients_main[1,1] * exp((coefficients_main[1,2]*x) * (p/coefficients_main[1,3]+p))
  return(predicted_ebu_rate)
}

SOA_coefficient_low <- function(x, p){
  predicted_ebu_rate = coefficients_low[1,1] * exp((coefficients_low[1,2]*x) * (p/coefficients_low[1,3]+p)) + model_error
  return(predicted_ebu_rate)
}

SOA_coefficient_high <- function(x, p){
  predicted_ebu_rate = coefficients_high[1,1] * exp((coefficients_high[1,2]*x) * (p/coefficients_high[1,3]+p)) - model_error
  return(predicted_ebu_rate)
}



LM_mean_ebu = glm(ch4_ebu ~ mean_sw_wm2, data = ebu_base_temp)

summary(LM_mean_ebu) # get MSE value

eval_models <- as.data.frame(cbind(ebu_base_temp$ch4_ebu, predict(LM_mean_ebu)))

eval_LM_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_mean <- rbind(LM_mean_ebu$coefficients)

coefficients_low <- cbind(LM_mean_ebu$coefficients[1]-summary(LM_mean_ebu)$coefficients[,2][1],
                          LM_mean_ebu$coefficients[2]-summary(LM_mean_ebu)$coefficients[,2][2])

coefficients_high <- cbind(LM_mean_ebu$coefficients[1]+summary(LM_mean_ebu)$coefficients[,2][1],
                          LM_mean_ebu$coefficients[2]+summary(LM_mean_ebu)$coefficients[,2][2])


LM_coefficient_mean <- function(x, p){
  predicted_ebu_rate = coefficients_mean[1,1] + p*coefficients_mean[1,2]
  return(predicted_ebu_rate)
}

LM_coefficient_low <- function(x, p){
  predicted_ebu_rate = coefficients_low[1,1] + p*coefficients_low[1,2]
  return(predicted_ebu_rate)
}

LM_coefficient_high <- function(x, p){
  predicted_ebu_rate = coefficients_high[1,1] + p*coefficients_high[1,2]
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

coefficients_main <- rbind(mean_ebu_MFOA_res$m$getPars())


coefficients_low <- cbind(mean_ebu_MFOA_res$m$getPars()[1]-summary(mean_ebu_MFOA_res)$coefficients[,2][1],
                          mean_ebu_MFOA_res$m$getPars()[2]-summary(mean_ebu_MFOA_res)$coefficients[,2][2])

coefficients_high <- cbind(mean_ebu_MFOA_res$m$getPars()[1]+summary(mean_ebu_MFOA_res)$coefficients[,2][1],
                           mean_ebu_MFOA_res$m$getPars()[2]+summary(mean_ebu_MFOA_res)$coefficients[,2][2])

FOA_coefficient_mean_res <- function(x){
  predicted_ebu_rate = (coefficients_main[1,1] * coefficients_main[1,2]^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_low_res <- function(x){
  predicted_ebu_rate = (coefficients_low[1,1] * coefficients_low[1,2]^(x-20))
  return(predicted_ebu_rate)
}

FOA_coefficient_high_res <- function(x){
  predicted_ebu_rate = (coefficients_high[1,1] * coefficients_high[1,2]^(x-20))
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

coefficients_main <- rbind(mean_ebu_SOA_temp_soil_threshold_res$m$getPars())


coefficients_low <- cbind(mean_ebu_SOA_temp_soil_threshold_res$m$getPars()[1]-summary(mean_ebu_SOA_temp_soil_threshold_res)$coefficients[,2][1],
                          mean_ebu_SOA_temp_soil_threshold_res$m$getPars()[2]-summary(mean_ebu_SOA_temp_soil_threshold_res)$coefficients[,2][2],
                          mean_ebu_SOA_temp_soil_threshold_res$m$getPars()[3]-summary(mean_ebu_SOA_temp_soil_threshold_res)$coefficients[,2][3])

coefficients_high <- cbind(mean_ebu_SOA_temp_soil_threshold_res$m$getPars()[1]+summary(mean_ebu_SOA_temp_soil_threshold_res)$coefficients[,2][1],
                           mean_ebu_SOA_temp_soil_threshold_res$m$getPars()[2]+summary(mean_ebu_SOA_temp_soil_threshold_res)$coefficients[,2][2],
                           mean_ebu_SOA_temp_soil_threshold_res$m$getPars()[3]+summary(mean_ebu_SOA_temp_soil_threshold_res)$coefficients[,2][3])

SOA_coefficient_mean_res <- function(x, p){
  predicted_ebu_rate = coefficients_main[1,1] * exp((coefficients_main[1,2]*x) * (p/coefficients_main[1,3]+p))
  return(predicted_ebu_rate)
}

SOA_coefficient_low_res <- function(x, p){
  predicted_ebu_rate = coefficients_low[1,1] * exp((coefficients_low[1,2]*x) * (p/coefficients_low[1,3]+p)) + model_error
  return(predicted_ebu_rate)
}

SOA_coefficient_high_res <- function(x, p){
  predicted_ebu_rate = coefficients_high[1,1] * exp((coefficients_high[1,2]*x) * (p/coefficients_high[1,3]+p)) - model_error
  return(predicted_ebu_rate)
}



LM_mean_ebu_res = glm(ch4_ebu ~ mean_sw_wm2, data = ebu_base_temp)

summary(LM_mean_ebu_res) # get MSE value

eval_models <- as.data.frame(cbind(ebu_base_temp$ch4_ebu, predict(LM_mean_ebu_res)))

eval_LM_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_mean <- rbind(LM_mean_ebu_res$coefficients)

coefficients_low <- cbind(LM_mean_ebu_res$coefficients[1]-summary(LM_mean_ebu_res)$coefficients[,2][1],
                          LM_mean_ebu_res$coefficients[2]-summary(LM_mean_ebu_res)$coefficients[,2][2])

coefficients_high <- cbind(LM_mean_ebu_res$coefficients[1]+summary(LM_mean_ebu_res)$coefficients[,2][1],
                           LM_mean_ebu_res$coefficients[2]+summary(LM_mean_ebu_res)$coefficients[,2][2])


LM_coefficient_mean_res <- function(x, p){
  predicted_ebu_rate = coefficients_mean[1,1] + p*coefficients_mean[1,2]
  return(predicted_ebu_rate)
}

LM_coefficient_low_res <- function(x, p){
  predicted_ebu_rate = coefficients_low[1,1] + p*coefficients_low[1,2]
  return(predicted_ebu_rate)
}

LM_coefficient_high_res <- function(x, p){
  predicted_ebu_rate = coefficients_high[1,1] + p*coefficients_high[1,2]
  return(predicted_ebu_rate)
}