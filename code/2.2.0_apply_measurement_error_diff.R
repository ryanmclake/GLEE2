library(minpack.lm)


### MONTHLY TIMESTEP LAKES ###
diff_base_temp <- error_diff %>%
  ungroup(.) %>%
  select(ch4_diff, sd_time, waterbody_type, sd_space, temp_for_model_C) %>%
  mutate(diff_time_low = ifelse(ch4_diff-sd_time < 0, 0.01, ch4_diff-sd_time)) %>%
  mutate(diff_time_high = ifelse(ch4_diff+sd_time < 0, 0.01, ch4_diff+sd_time)) %>%
  mutate(diff_space_low = ifelse(ch4_diff-sd_space < 0, 0.01, ch4_diff-sd_space)) %>%
  mutate(diff_space_high = ifelse(ch4_diff+sd_space < 0, 0.01, ch4_diff+sd_space)) %>%
  na.omit(.)%>%
  filter(waterbody_type == "lake")

# First-order Arhennius
mean_diff_FOA = nlsLM(ch4_diff ~ A * exp(a * temp_for_model_C),
                        start = list(A = 0.023, a = 0.124),
                        data = diff_base_temp,
                        control = nls.lm.control(maxiter=1000))


time_low_diff_FOA = nlsLM(diff_time_low ~ A * exp(a * temp_for_model_C),
                      start = list(A = 0.023, a = 0.124),
                      data = diff_base_temp,
                      control = nls.lm.control(maxiter=1000))

time_high_diff_FOA = nlsLM(diff_time_high ~ A * exp(a * temp_for_model_C),
                      start = list(A = 0.023, a = 0.124),
                      data = diff_base_temp,
                      control = nls.lm.control(maxiter=1000))

space_low_diff_FOA = nlsLM(diff_space_low ~ A * exp(a * temp_for_model_C),
                      start = list(A = 0.023, a = 0.124),
                      data = diff_base_temp,
                      control = nls.lm.control(maxiter=1000))

space_high_diff_FOA = nlsLM(diff_space_high ~ A * exp(a * temp_for_model_C),
                      start = list(A = 0.023, a = 0.124),
                      data = diff_base_temp,
                      control = nls.lm.control(maxiter=1000))


summary(mean_diff_FOA) # get MSE value
summary(time_low_diff_FOA) # get MSE value
summary(time_high_diff_FOA) # get MSE value
summary(space_low_diff_FOA) # get MSE value
summary(space_high_diff_FOA) # get MSE value

eval_models <- as.data.frame(cbind(diff_base_temp$ch4_diff, predict(mean_diff_FOA), predict(time_low_diff_FOA),predict(time_high_diff_FOA),
                                   predict(space_low_diff_FOA), predict(space_high_diff_FOA)))

eval_FOA_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1),
            NSE_time_low = NSE(V3, V1),
            NSE_time_high = NSE(V4, V1),
            NSE_space_low = NSE(V5, V1),
            NSE_space_high = NSE(V6, V1))

coefficients_a <- rbind(mean_diff_FOA$m$getPars(),time_low_diff_FOA$m$getPars(),
                        time_high_diff_FOA$m$getPars(),space_low_diff_FOA$m$getPars(),
                        space_high_diff_FOA$m$getPars())

FOA_mean <- function(x){
  predicted_diff_rate = 7.1983292 * exp(0.06377818*x)
  return(predicted_diff_rate)
}

FOA_time_low <-  function(x){
  predicted_diff_rate = 2.4784471 * exp(0.07316942*x)
  return(predicted_diff_rate)
}

FOA_time_high <-  function(x){
  predicted_diff_rate = 63.8196323 * exp(0.03854356*x)
  return(predicted_diff_rate)
}

FOA_space_low <-  function(x){
  predicted_diff_rate = 0.1420812 * exp(0.16359689*x)
  return(predicted_diff_rate)
}

FOA_space_high <-  function(x){
  predicted_diff_rate = 89.4296262 * exp(0.03529463*x)
  return(predicted_diff_rate)
}


### LAKE AREA ALL DATA ###

diff_base_temp <- error_diff %>%
  ungroup(.) %>%
  select(ch4_diff, sd_time, waterbody_type, sd_space, surf_area_k) %>%
  mutate(diff_time_low = ifelse(ch4_diff-sd_time < 0, 0.01, ch4_diff-sd_time)) %>%
  mutate(diff_time_high = ifelse(ch4_diff+sd_time < 0, 0.01, ch4_diff+sd_time)) %>%
  mutate(diff_space_low = ifelse(ch4_diff-sd_space < 0, 0.01, ch4_diff-sd_space)) %>%
  mutate(diff_space_high = ifelse(ch4_diff+sd_space < 0, 0.01, ch4_diff+sd_space)) %>%
  na.omit(.)%>%
  filter(waterbody_type == "lake")


LM_mean_diff = glm(ch4_diff ~ surf_area_k, data = diff_base_temp)
LM_time_low_diff = glm(diff_time_low ~ surf_area_k, data = diff_base_temp)
LM_time_high_diff = glm(diff_time_high ~ surf_area_k, data = diff_base_temp)
LM_space_low_diff = glm(diff_space_low ~ surf_area_k, data = diff_base_temp)
LM_space_high_diff = glm(diff_space_high ~ surf_area_k, data = diff_base_temp)

summary(LM_mean_diff) # get MSE value
summary(LM_time_low_diff) # get MSE value
summary(LM_time_high_diff) # get MSE value
summary(LM_space_low_diff) # get MSE value
summary(LM_space_high_diff) # get MSE value

eval_models <- as.data.frame(cbind(diff_base_temp$ch4_diff, predict(LM_mean_diff), predict(LM_time_low_diff),predict(LM_time_high_diff),
                                   predict(LM_space_low_diff), predict(LM_space_high_diff)))


eval_LM_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1),
            NSE_time_low = NSE(V3, V1),
            NSE_time_high = NSE(V4, V1),
            NSE_space_low = NSE(V5, V1),
            NSE_space_high = NSE(V6, V1))


coefficients_b <- rbind(LM_mean_diff$coefficients,LM_time_low_diff$coefficients,
                        LM_time_high_diff$coefficients,LM_space_low_diff$coefficients,
                        LM_space_high_diff$coefficients)


LM_mean <- function(x, p){
  predicted_diff_rate = coefficients_b[1,1] + p*coefficients_b[1,2]
  return(predicted_diff_rate)
}

LM_time_low <- function(x, p){
  predicted_diff_rate = coefficients_b[2,1] + p*coefficients_b[2,2]
  return(predicted_diff_rate)
}

LM_time_high <- function(x, p){
  predicted_diff_rate = coefficients_b[3,1] + p*coefficients_b[3,2]
  return(predicted_diff_rate)
}

LM_space_low <- function(x, p){
  predicted_diff_rate = coefficients_b[4,1] + p*coefficients_b[4,2]
  return(predicted_diff_rate)
}

LM_space_high <- function(x, p){
  predicted_diff_rate = coefficients_b[5,1] + p*coefficients_b[5,2]
  return(predicted_diff_rate)
}





### RESERVOIRS ###


### MONTHLY TIMESTEP LAKES ###
diff_base_temp <- error_diff %>%
  ungroup(.) %>%
  select(ch4_diff, sd_time, waterbody_type, sd_space, temp_for_model_C) %>%
  mutate(diff_time_low = ifelse(ch4_diff-sd_time < 0, 0.01, ch4_diff-sd_time)) %>%
  mutate(diff_time_high = ifelse(ch4_diff+sd_time < 0, 0.01, ch4_diff+sd_time)) %>%
  mutate(diff_space_low = ifelse(ch4_diff-sd_space < 0, 0.01, ch4_diff-sd_space)) %>%
  mutate(diff_space_high = ifelse(ch4_diff+sd_space < 0, 0.01, ch4_diff+sd_space)) %>%
  na.omit(.) %>%
  filter(waterbody_type == "reservoir")

# First-order Arhennius
mean_diff_FOA_res = nlsLM(ch4_diff ~ A * exp(a * temp_for_model_C),
                      start = list(A = 0.023, a = 0.124),
                      data = diff_base_temp,
                      control = nls.lm.control(maxiter=1000))


time_low_diff_FOA_res = nlsLM(diff_time_low ~ A * exp(a * temp_for_model_C),
                          start = list(A = 0.023, a = 0.124),
                          data = diff_base_temp,
                          control = nls.lm.control(maxiter=1000))

time_high_diff_FOA_res = nlsLM(diff_time_high ~ A * exp(a * temp_for_model_C),
                           start = list(A = 0.023, a = 0.124),
                           data = diff_base_temp,
                           control = nls.lm.control(maxiter=1000))

space_low_diff_FOA_res = nlsLM(diff_space_low ~ A * exp(a * temp_for_model_C),
                           start = list(A = 0.023, a = 0.124),
                           data = diff_base_temp,
                           control = nls.lm.control(maxiter=1000))

space_high_diff_FOA_res = nlsLM(diff_space_high ~ A * exp(a * temp_for_model_C),
                            start = list(A = 0.023, a = 0.124),
                            data = diff_base_temp,
                            control = nls.lm.control(maxiter=1000))


summary(monthly_diff_FOA_res) # get MSE value
summary(time_low_diff_FOA_res) # get MSE value
summary(time_high_diff_FOA_res) # get MSE value
summary(space_low_diff_FOA_res) # get MSE value
summary(space_high_diff_FOA_res) # get MSE value

eval_models <- as.data.frame(cbind(diff_base_temp$ch4_diff, predict(mean_diff_FOA_res), predict(time_low_diff_FOA_res),predict(time_high_diff_FOA_res),
                                   predict(space_low_diff_FOA_res), predict(space_high_diff_FOA_res)))

eval_FOA_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1),
            NSE_time_low = NSE(V3, V1),
            NSE_time_high = NSE(V4, V1),
            NSE_space_low = NSE(V5, V1),
            NSE_space_high = NSE(V6, V1))

coefficients_a <- rbind(mean_diff_FOA_res$m$getPars(),time_low_diff_FOA_res$m$getPars(),
                        time_high_diff_FOA_res$m$getPars(),space_low_diff_FOA_res$m$getPars(),
                        space_high_diff_FOA_res$m$getPars())

FOA_mean_res <- function(x){
  predicted_diff_rate = 14.507176 * exp(0.02903751*x)
  return(predicted_diff_rate)
}

FOA_time_low_res <-  function(x){
  predicted_diff_rate = 7.429610 * exp(0.01730997*x)
  return(predicted_diff_rate)
}

FOA_time_high_res <-  function(x){
  predicted_diff_rate = 394.371565 * exp(-0.05280574*x)
  return(predicted_diff_rate)
}

FOA_space_low_res <-  function(x){
  predicted_diff_rate = 6.151718 * exp(0.02983359*x)
  return(predicted_diff_rate)
}

FOA_space_high_res <-  function(x){
  predicted_diff_rate = 562.518005 * exp(-0.06166063*x)
  return(predicted_diff_rate)
}


### RES AREA ALL DATA ###

diff_base_temp <- error_diff %>%
  ungroup(.) %>%
  select(ch4_diff, sd_time, waterbody_type, sd_space, surf_area_k) %>%
  mutate(diff_time_low = ifelse(ch4_diff-sd_time < 0, 0.01, ch4_diff-sd_time)) %>%
  mutate(diff_time_high = ifelse(ch4_diff+sd_time < 0, 0.01, ch4_diff+sd_time)) %>%
  mutate(diff_space_low = ifelse(ch4_diff-sd_space < 0, 0.01, ch4_diff-sd_space)) %>%
  mutate(diff_space_high = ifelse(ch4_diff+sd_space < 0, 0.01, ch4_diff+sd_space)) %>%
  na.omit(.)%>%
  filter(waterbody_type == "reservoir")


LM_mean_diff_res = glm(ch4_diff ~ surf_area_k, data = diff_base_temp)
LM_time_low_diff_res = glm(diff_time_low ~ surf_area_k, data = diff_base_temp)
LM_time_high_diff_res = glm(diff_time_high ~ surf_area_k, data = diff_base_temp)
LM_space_low_diff_res = glm(diff_space_low ~ surf_area_k, data = diff_base_temp)
LM_space_high_diff_res = glm(diff_space_high ~ surf_area_k, data = diff_base_temp)

summary(LM_mean_diff_res) # get MSE value
summary(LM_time_low_diff_res) # get MSE value
summary(LM_time_high_diff_res) # get MSE value
summary(LM_space_low_diff_res) # get MSE value
summary(LM_space_high_diff_res) # get MSE value

eval_models <- as.data.frame(cbind(diff_base_temp$ch4_diff, predict(LM_mean_diff_res), predict(LM_time_low_diff_res),predict(LM_time_high_diff_res),
                                   predict(LM_space_low_diff_res), predict(LM_space_high_diff_res)))


eval_LM_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1),
            NSE_time_low = NSE(V3, V1),
            NSE_time_high = NSE(V4, V1),
            NSE_space_low = NSE(V5, V1),
            NSE_space_high = NSE(V6, V1))


coefficients_b <- rbind(LM_mean_diff_res$coefficients,LM_time_low_diff_res$coefficients,
                        LM_time_high_diff_res$coefficients,LM_space_low_diff_res$coefficients,
                        LM_space_high_diff_res$coefficients)


LM_mean_res <- function(x, p){
  predicted_diff_rate = coefficients_b[1,1] + p*coefficients_b[1,2]
  return(predicted_diff_rate)
}

LM_time_low_res <- function(x, p){
  predicted_diff_rate = coefficients_b[2,1] + p*coefficients_b[2,2]
  return(predicted_diff_rate)
}

LM_time_high_res <- function(x, p){
  predicted_diff_rate = coefficients_b[3,1] + p*coefficients_b[3,2]
  return(predicted_diff_rate)
}

LM_space_low_res <- function(x, p){
  predicted_diff_rate = coefficients_b[4,1] + p*coefficients_b[4,2]
  return(predicted_diff_rate)
}

LM_space_high_res <- function(x, p){
  predicted_diff_rate = coefficients_b[5,1] + p*coefficients_b[5,2]
  return(predicted_diff_rate)
}

























# 
# 
# diff_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
#                                   values = c(
#                                     monthly_FOA_function(-15:35),
#                                     natchimuthu_FOA_function(-15:35)),
#                                   model = rep(c(
#                                     "First-Order (F-O) Arhennius",
#                                     "Natchimuthu's F-O Arhennius"), each = 51))
# 
# ggplot(diff_base_temp, aes(x = temp_for_model_C, y = ch4_diff))+
#   geom_point()+
#   geom_line(data = diff_function_output, aes(x, values, group = model, color = model), lwd = 2)+
#   scale_color_viridis(discrete = T, option = "C")+
#   coord_cartesian(ylim=c(0, 2000))+
#   theme_classic()
# 
# 
# 
# ### YEARLY TIMESTEP ###
# 
# diff_base_yearly <- base %>% select(year, ch4_diff, waterbody_id, temp_for_model_K) %>%
#   na.omit(.)%>%
#   mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
#   group_by(waterbody_id, year)%>%
#   summarize_all(funs(mean), na.rm = T)
# 
# # First-order Arhennius
# yearly_diff_FOA = nlsLM(ch4_diff ~ A * exp(a * temp_for_model_C),
#                        start = list(A = 0.12, a = 0.035),
#                        data = diff_base_yearly,
#                        control = nls.lm.control(maxiter=1000))
# summary(yearly_diff_FOA) # get MSE value
# 
# 
# yearly_FOA_function <- function(x) {18.05549 * exp(0.03913*x)}
# natchimuthu_FOA_function <- function(x) {0.023 * exp(0.124*x)}
# 
# diff_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
#                                    values = c(
#                                      yearly_FOA_function(-15:35),
#                                      natchimuthu_FOA_function(-15:35)),
#                                    model = rep(c(
#                                      "First-Order (F-O) Arhennius",
#                                      "Natchimuthu's F-O Arhennius"), each = 51))
# 
# ggplot(diff_base_yearly, aes(x = temp_for_model_C, y = ch4_diff))+
#   geom_point()+
#   geom_line(data = diff_function_output, aes(x, values, group = model, color = model), lwd = 2)+
#   scale_color_viridis(discrete = T, option = "C")+
#   coord_cartesian(ylim=c(0, 3000))+
#   theme_classic()
# 
# ### SITE MEAN ###
# 
# diff_base_site <- base %>% select(ch4_diff, waterbody_id, temp_for_model_K) %>%
#   na.omit(.)%>%
#   mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
#   group_by(waterbody_id)%>%
#   summarize_all(funs(mean), na.rm = T)
# 
# # First-order Arhennius
# site_diff_FOA = nlsLM(ch4_diff ~ A * exp(a * temp_for_model_C),
#                      start = list(A = 0.12, a = 0.035),
#                      data = diff_base_site,
#                      control = nls.lm.control(maxiter=1000))
# summary(site_diff_FOA) # get MSE value
# 
# site_FOA_function <- function(x) {19.37369 * exp(0.03704*x)}
# natchimuthu_FOA_function <- function(x) {0.023 * exp(0.124*x)}
# 
# diff_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
#                                    values = c(
#                                      yearly_FOA_function(-15:35),
#                                      natchimuthu_FOA_function(-15:35)),
#                                    model = rep(c(
#                                      "First-Order (F-O) Arhennius",
#                                      "Natchimuthu's F-O Arhennius"), each = 51))
# 
# ggplot(diff_base_site, aes(x = temp_for_model_C, y = ch4_diff))+
#   geom_point()+
#   geom_line(data = diff_function_output, aes(x, values, group = model, color = model), lwd = 2)+
#   scale_color_viridis(discrete = T, option = "C")+
#   coord_cartesian(ylim=c(0, 3000))+
#   theme_classic()
# 
# 
# ### LAKE AREA MODELS ###
# 
# 
# 
# 
# 
