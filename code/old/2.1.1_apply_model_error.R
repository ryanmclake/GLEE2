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

coefficients_a <- rbind(mean_ebu_MFOA$m$getPars())
model_error <- sqrt(deviance(mean_ebu_MFOA))

FOA_mean <- function(x){
  predicted_ebu_rate = (coefficients_a[1,1] * coefficients_a[1,2]^(x-20))
  return(predicted_ebu_rate)
}

FOA_mean_error_low <- function(x){
  predicted_ebu_rate = (coefficients_a[1,1] * coefficients_a[1,2]^(x-20)) - model_error
  return(predicted_ebu_rate)
}

FOA_mean_error_high <- function(x){
  predicted_ebu_rate = (coefficients_a[1,1] * coefficients_a[1,2]^(x-20)) + model_error
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

coefficients_b <- rbind(mean_ebu_SOA_temp_soil_threshold$m$getPars())
model_error <- sqrt(deviance(mean_ebu_SOA_temp_soil_threshold))


SOA_mean <- function(x, p){
  predicted_ebu_rate = coefficients_b[1,1] * exp((coefficients_b[1,2]*x) * (p/coefficients_b[1,3]+p))
  return(predicted_ebu_rate)
}

SOA_error_low <- function(x, p){
  predicted_ebu_rate = coefficients_b[1,1] * exp((coefficients_b[1,2]*x) * (p/coefficients_b[1,3]+p)) + model_error
  return(predicted_ebu_rate)
}

SOA_error_high <- function(x, p){
  predicted_ebu_rate = coefficients_b[1,1] * exp((coefficients_b[1,2]*x) * (p/coefficients_b[1,3]+p)) - model_error
  return(predicted_ebu_rate)
}



LM_mean_ebu = glm(ch4_ebu ~ mean_sw_wm2, data = ebu_base_temp)

summary(LM_mean_ebu) # get MSE value

eval_models <- as.data.frame(cbind(ebu_base_temp$ch4_ebu, predict(LM_mean_ebu)))

eval_LM_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_b <- rbind(LM_mean_ebu$coefficients)
model_error <- sqrt(deviance(LM_mean_ebu))


LM_mean <- function(x, p){
  predicted_ebu_rate = coefficients_b[1,1] + p*coefficients_b[1,2]
  return(predicted_ebu_rate)
}

LM_error_low <- function(x, p){
  predicted_ebu_rate = coefficients_b[1,1] + p*coefficients_b[1,2] - model_error
  return(predicted_ebu_rate)
}

LM_error_high <- function(x, p){
  predicted_ebu_rate = coefficients_b[1,1] + p*coefficients_b[1,2] + model_error
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

coefficients_a <- rbind(mean_ebu_MFOA_res$m$getPars())
model_error <- sqrt(deviance(mean_ebu_MFOA_res))

FOA_mean <- function(x){
  predicted_ebu_rate = (coefficients_a[1,1] * coefficients_a[1,2]^(x-20))
  return(predicted_ebu_rate)
}

FOA_mean_error_low <- function(x){
  predicted_ebu_rate = (coefficients_a[1,1] * coefficients_a[1,2]^(x-20)) - model_error
  return(predicted_ebu_rate)
}

FOA_mean_error_high <- function(x){
  predicted_ebu_rate = (coefficients_a[1,1] * coefficients_a[1,2]^(x-20)) + model_error
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

coefficients_b <- rbind(mean_ebu_SOA_temp_soil_threshold_res$m$getPars())
model_error <- sqrt(deviance(mean_ebu_SOA_temp_soil_threshold_res))


SOA_mean <- function(x, p){
  predicted_ebu_rate = coefficients_b[1,1] * exp((coefficients_b[1,2]*x) * (p/coefficients_b[1,3]+p))
  return(predicted_ebu_rate)
}

SOA_error_low <- function(x, p){
  predicted_ebu_rate = coefficients_b[1,1] * exp((coefficients_b[1,2]*x) * (p/coefficients_b[1,3]+p)) + model_error
  return(predicted_ebu_rate)
}

SOA_error_high <- function(x, p){
  predicted_ebu_rate = coefficients_b[1,1] * exp((coefficients_b[1,2]*x) * (p/coefficients_b[1,3]+p)) - model_error
  return(predicted_ebu_rate)
}


LM_mean_ebu_res = glm(ch4_ebu ~ mean_sw_wm2, data = ebu_base_temp)

summary(LM_mean_ebu_res) # get MSE value

eval_models <- as.data.frame(cbind(ebu_base_temp$ch4_ebu, predict(LM_mean_ebu_res)))

eval_LM_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_b <- rbind(LM_mean_ebu_res$coefficients)
model_error <- sqrt(deviance(LM_mean_ebu_res))


LM_mean <- function(x, p){
  predicted_ebu_rate = coefficients_b[1,1] + p*coefficients_b[1,2]
  return(predicted_ebu_rate)
}

LM_error_low <- function(x, p){
  predicted_ebu_rate = coefficients_b[1,1] + p*coefficients_b[1,2] - model_error
  return(predicted_ebu_rate)
}

LM_error_high <- function(x, p){
  predicted_ebu_rate = coefficients_b[1,1] + p*coefficients_b[1,2] + model_error
  return(predicted_ebu_rate)
}












































# 
# # # fake some reference data
# # ref<-ebu_base_temp$ch4_ebu
# # # add a little noise
# # model1<-predict(monthly_ebu_MFOA)
# # model2<-predict(monthly_ebu_SOA_temp_threshold)
# # # add more noise
# # model2<-ref+rnorm(30)
# # # display the diagram with the better model
# # oldpar<-taylor.diagram(ref,model1)
# # # now add the worse model
# # taylor.diagram(ref,model1,add=TRUE,col="blue")
# # # get approximate legend position
# # lpos<-1.5*sd(ref)
# # # add a legend
# # legend(lpos,lpos,legend=c("Better","Worse"),pch=19,col=c("red","blue"))
# # # now restore par values
# # par(oldpar)
# # # show the "all correlation" display
# # taylor.diagram(ref,model1,pos.cor=T)
# # taylor.diagram(ref,model2,add=TRUE,col="blue")
# 
# 
# # 
# # ebu_base_g_res <- ebu_base_models %>%
# #   select(ch4_ebu, mean_sw_wm2) %>%
# #   na.omit(.) %>%
# #   mutate(littoral_area = 40)
# # 
# # # G-Res Models
# # ebullition_base_model <- function(mean_sw_wm2, littoral_area){
# #   10^(-0.98574 + 1.0075 * log10(littoral_area) + 0.04928 * ((mean_sw_wm2)/30.4))
# # }
# # 
# # ebullition_base_model <- Vectorize(ebullition_base_model)
# # 
# # g_res <- ebu_base_g_res %>%
# #   mutate(flux = ebullition_base_model(mean_sw_wm2 = mean_sw_wm2, littoral_area = littoral_area))
# # #Where --> 
# # # Parameters
# # # B1 = -0.98574
# # # B2 = 1.0075
# # # B3 = 0.04928
# # 
# # # Inputs - covariates
# # # per_littoral_area = littoral fraction above 3 meters
# # # I_cumu = Cumulative global Hiorizontal radiance 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# roundUp <- function(x,to=3){
#   to*(x%/%to + as.logical(x%%to))
# }
# 
# ebu_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
#                                   values = c(
#                                    monthly_MFOA_function(-15:35),
#                                    aben_MFOA_function(-15:35),
#                                    monthly_ebu_SOA_function(-15:35)),
#                                   model = rep(c(
#                                     "First-Order Arhennius",
#                                     "Aben's Arhennius",
#                                     "Second-Order Exponential"), each = 51))
# 
# error_ebu2 <- error_ebu %>%
#   filter(waterbody_type == "reservoir") %>%
#   mutate(temp_rnd = roundUp(temp_for_model_K)) %>%
#   group_by(temp_rnd) %>%
#   summarise(mean_ebu = mean(ch4_ebu),
#             sd_ebu = sd(ch4_ebu))
# 
# ggplot(error_ebu2, aes(x = temp_rnd-273.15, y = mean_ebu))+
#   geom_point()+
#   geom_errorbar(aes(x = temp_rnd-273.15, y = mean_ebu, ymin = mean_ebu-sd_ebu, ymax = mean_ebu+sd_ebu))+
#   geom_line(data = ebu_function_output, aes(x, values, group = model, color = model), lwd = 2, alpha = 0.4)+
#   scale_color_viridis(discrete = T, option = "C")+
#   coord_cartesian(ylim=c(0, 1700))+
#   theme_classic()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ### YEARLY TIMESTEP ###
# 
# ebu_base_yearly <- base %>% select(year, ch4_ebu, waterbody_id, temp_for_model_K) %>%
#   na.omit(.)%>%
#   mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
#   group_by(waterbody_id, year)%>%
#   summarize_all(funs(mean), na.rm = T)
# 
# # First-order Arhennius
# yearly_ebu_FOA = nlsLM(ch4_ebu ~ A * exp(a * temp_for_model_C),
#                         start = list(A = 0.12, a = 0.035),
#                         data = ebu_base_yearly,
#                         control = nls.lm.control(maxiter=1000))
# summary(yearly_ebu_FOA) # get MSE value
# 
# # Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
# yearly_ebu_MFOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
#                          start = list(A = 100, a = 1.1),
#                          data = ebu_base_yearly,
#                          control = nls.lm.control(maxiter=1000))
# summary(yearly_ebu_MFOA) # get MSE value
# 
# # Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
# yearly_ebu_MFOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
#                          start = list(A = 100, a = 1.1),
#                          data = ebu_base_yearly,
#                          control = nls.lm.control(maxiter=1000))
# summary(yearly_ebu_MFOA) # get MSE value
# 
# yearly_FOA_function <- function(x) {13.71050 * exp(0.08813*x)}
# yearly_MFOA_function <- function(x) {79.89646 * 1.09213 ^ (x-20)}
# aben_MFOA_function <- function(x) {100 * 1.1 ^ (x-20)}
# 
# ebu_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
#                                   values = c(
#                                     yearly_FOA_function(-15:35),
#                                     yearly_MFOA_function(-15:35),
#                                     aben_MFOA_function(-15:35)),
#                                   model = rep(c(
#                                     "First-Order (F-O) Arhennius",
#                                     "Modified F-O Arhennius",
#                                     "Aben's Modified F-O Arhennius"), each = 51))
# 
# ggplot(ebu_base_yearly, aes(x = temp_for_model_C, y = ch4_ebu))+
#   geom_point()+
#   geom_line(data = ebu_function_output, aes(x, values, group = model, color = model), lwd = 2)+
#   scale_color_viridis(discrete = T, option = "C")+
#   coord_cartesian(ylim=c(0, 3000))+
#   theme_classic()
# 
# ### SITE MEAN ###
# 
# ebu_base_site <- base %>% select(ch4_ebu, waterbody_id, temp_for_model_K) %>%
#   na.omit(.)%>%
#   mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
#   group_by(waterbody_id)%>%
#   summarize_all(funs(mean), na.rm = T)
# 
# # First-order Arhennius
# site_ebu_FOA = nlsLM(ch4_ebu ~ A * exp(a * temp_for_model_C),
#                        start = list(A = 0.12, a = 0.035),
#                        data = ebu_base_site,
#                        control = nls.lm.control(maxiter=1000))
# summary(site_ebu_FOA) # get MSE value
# 
# # Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
# site_ebu_MFOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
#                         start = list(A = 100, a = 1.1),
#                         data = ebu_base_site,
#                         control = nls.lm.control(maxiter=1000))
# summary(site_ebu_MFOA) # get MSE value
# 
# site_FOA_function <- function(x) {19.18062 * exp(0.06981*x)}
# site_MFOA_function <- function(x) {77.48887 * 1.07230 ^ (x-20)}
# aben_MFOA_function <- function(x) {100 * 1.1 ^ (x-20)}
# 
# ebu_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
#                                   values = c(
#                                     site_FOA_function(-15:35),
#                                     site_MFOA_function(-15:35),
#                                     aben_MFOA_function(-15:35)),
#                                   model = rep(c(
#                                     "First-Order (F-O) Arhennius",
#                                     "Modified F-O Arhennius",
#                                     "Aben's Modified F-O Arhennius"), each = 51))
# 
# ggplot(ebu_base_site, aes(x = temp_for_model_C, y = ch4_ebu))+
#   geom_point()+
#   geom_line(data = ebu_function_output, aes(x, values, group = model, color = model), lwd = 2)+
#   scale_color_viridis(discrete = T, option = "C")+
#   coord_cartesian(ylim=c(0, 3000))+
#   theme_classic()
# 
# 
# ### LAKE AREA MODELS ###
# 
# write_csv(base2, "./data/GLEE_data_with_GLCP_HWSD_link.csv")
# 
# base2 <- read_csv("./data/GLEE_data_with_GLCP_HWSD_link.csv")
# 
# ebu_base_temp_soil <- base2 %>% select(ch4_ebu, temp_for_model_K, waterbody_id, month, T_OC) %>%
#   na.omit(.)%>%
#   mutate(temp_for_model_C = temp_for_model_K-273.15)
# 
# 
# monthly_ebu_SOA_soil_threshold = nlsLM(ch4_ebu ~ (A * exp(a * temp_for_model_C)) * (-exp(b * T_OC)),
#                                        start = list(A = 0.12, a = 0.035, b = 0.12),
#                                        data = ebu_base_temp_soil,
#                                        control = nls.lm.control(maxiter=1024))
# summary(monthly_ebu_SOA_soil_threshold) # get MSE value
# 
# dat <- as.data.frame(cbind(ebu_base_temp_soil$ch4_ebu, predict(monthly_ebu_SOA_soil_threshold)))
# second_order_area_limiter_NSE <- NSE(dat$V2, dat$V1)
# 
# 
# ### LAKE AREA ALL DATA ###
# 
# ebu_base_area <- base %>% select(ch4_ebu, waterbody_id, month, surf_area_k) %>%
#   na.omit(.) %>%
#   mutate(ch4_ebu_scaled = scale(ch4_ebu),
#          surf_area_k_scaled = scale(surf_area_k))
# 
# area_monthly_model <- lm(ch4_ebu_scaled ~ surf_area_k_scaled, data = ebu_base_area)
# summary(area_monthly_model)
# 
# 
# 
# 
# ebu_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
#                                   values = c(
#                                     monthly_MFOA_function(-15:35),
#                                     yearly_MFOA_function(-15:35),
#                                     site_MFOA_function(-15:35),
#                                     aben_MFOA_function(-15:35)),
#                                   model = rep(c(
#                                     "Monthly First-Order Arhennius",
#                                     "Annual First-Order Arhennius",
#                                     "Site-Level First-Order Arhennius",
#                                     "Aben's Modified F-O Arhennius"), each = 51))
# 
# ggplot(ebu_base_temp, aes(x = temp_for_model_C, y = ch4_ebu))+
#   geom_point()+
#   geom_line(data = ebu_function_output, aes(x, values, group = model, color = model), lwd = 2)+
#   scale_color_viridis(discrete = T, option = "C")+
#   coord_cartesian(ylim=c(0, 3000))+
#   theme_classic()
# 
# 
# ### Lake and reservoir partition ###
# 
# ebu_base <- base %>% select(ch4_ebu, temp_for_model_K, waterbody_type, month) %>%
#   na.omit(.)%>%
#   mutate(temp_for_model_C = temp_for_model_K-273.15,
#          waterbody_type = ifelse(waterbody_type == "pond", "lake", waterbody_type)) 
# 
# ebu_base_reservoir <- base %>% select(ch4_ebu, temp_for_model_K, waterbody_type, month) %>%
#   na.omit(.)%>%
#   mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
#   filter(waterbody_type == "reservoir")
# 
# ebu_base_lake <- base %>% select(ch4_ebu, temp_for_model_K, waterbody_type, month) %>%
#   na.omit(.)%>%
#   mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
#   filter(waterbody_type == "lake")
# 
# # Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
# reservoir_ebu_FOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
#                          start = list(A = 100, a = 1.1),
#                          data = ebu_base_reservoir,
#                          control = nls.lm.control(maxiter=1000))
# summary(reservoir_ebu_FOA) # get MSE value
# 
# 
# 
# # Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
# lake_ebu_FOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
#                          start = list(A = 100, a = 1.1),
#                          data = ebu_base_lake,
#                          control = nls.lm.control(maxiter=1000))
# summary(lake_ebu_FOA) # get MSE value
# 
# 
# lake_FOA_function <- function(x) {77.07198 * 1.04613 ^ (x-20)}
# reservoir_MFOA_function <- function(x) {10.3765 * 1.6844 ^ (x-20)}
# aben_MFOA_function <- function(x) {100 * 1.1 ^ (x-20)}
# 
# ebu_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
#                                   values = c(
#                                     lake_FOA_function(-15:35),
#                                     reservoir_MFOA_function(-15:35),
#                                     aben_MFOA_function(-15:35)),
#                                   model = rep(c(
#                                     "First-Order Arhennius (Lake)",
#                                     "First-Order Arhennius (Reservoir)",
#                                     "Aben's Arhennius"), each = 51))
# 
# ggplot(ebu_base, aes(x = temp_for_model_C, y = ch4_ebu))+
#   geom_point(aes(shape = waterbody_type), size = 3)+
#   geom_line(data = ebu_function_output, aes(x, values, group = model, color = model), lwd = 2)+
#   scale_color_viridis(discrete = T, option = "C")+
#   coord_cartesian(ylim=c(0, 3000))+
#   theme_classic()
# 
# 
# 
# function1 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 0.001^2))} 
# function2 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 0.01^2))} 
# function3 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 0.1^2))} 
# function4 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 1^2))} 
# function5 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 10^2))} 
# function6 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 100^2))} 
# function7 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 1000^2))} 
# function8 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 2000^2))} 
# function9 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 3000^2))} 
# function10 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 4000^2))} 
# function11 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 5000^2))} 
# 
# 
# 
# 
# function_output <- data.frame(x = -15:35,            # Create data for ggplot2
#                               values = c(function1(-15:35),
#                                          function2(-15:35),
#                                          function3(-15:35),
#                                          function4(-15:35),
#                                          function5(-15:35),
#                                          function6(-15:35),
#                                          function7(-15:35),
#                                          function8(-15:35),
#                                          function9(-15:35),
#                                          function10(-15:35),
#                                          function11(-15:35)),
#                               surface_area = rep(c(0.001,
#                                                 0.01,
#                                                 0.1,
#                                                 1,
#                                                 10,
#                                                 100,
#                                                 1000,
#                                                 2000,
#                                                 3000,
#                                                 4000,
#                                                 5000), each = 51)) %>%
#   rename(`Lake Surface Area (km2)` = surface_area)
# 
# 
# 
# function_plot <- ggplot(function_output,
#                         aes(x, values, group = `Lake Surface Area (km2)`, col = `Lake Surface Area (km2)`)) +
#   geom_line()+
#   scale_color_viridis(discrete = F, option = "C")+
#   geom_point(data = ebu_base_temp_area, aes(x = temp_for_model_C, y = ch4_ebu), inherit.aes = F)+
#   theme_light()+
#   ylim(c(0,3000))+
#   labs(x = "Water Temperature (C)", y = "Ebullition Rate")
# 
# ggsave(function_plot, path = ".",
#        filename = "./figures/temp_do_elevation_concept_plot.jpg",
#        width = 8, height = 6, device='jpg', dpi=1000)
