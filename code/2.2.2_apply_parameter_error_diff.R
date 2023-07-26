library(minpack.lm)
library(hydroGOF)

set.seed(1098)

#How does observation, coefficient, and model variability contribute to total global emissions uncertainty?
### EXPLORE ERROR IN THE OBSERVATIONS FOR THE Q10 TEMP MODELS ###
# Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)

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


summary(mean_diff_FOA) # get MSE value

eval_models <- as.data.frame(cbind(diff_base_temp$ch4_diff, predict(mean_diff_FOA)))

eval_FOA_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_main <- rbind(mean_diff_FOA$m$getPars())


coefficients_low <- cbind(mean_diff_FOA$m$getPars()[1]-summary(mean_diff_FOA)$coefficients[,2][1],
                          mean_diff_FOA$m$getPars()[2]-summary(mean_diff_FOA)$coefficients[,2][2])

coefficients_high <- cbind(mean_diff_FOA$m$getPars()[1]+summary(mean_diff_FOA)$coefficients[,2][1],
                           mean_diff_FOA$m$getPars()[2]+summary(mean_diff_FOA)$coefficients[,2][2])

FOA_mean <- function(x){
  predicted_diff_rate = (coefficients_main[1,1] * exp(coefficients_main[1,2]*x))
  return(predicted_diff_rate)
}

FOA_coefficient_low_lake <- function(x){
  predicted_diff_rate = (3.907062 * exp(0.04560661*x))
  return(predicted_diff_rate)
}

FOA_coefficient_high_lake <- function(x){
  predicted_diff_rate = (10.4896 * exp(0.08194976*x))
  return(predicted_diff_rate)
}


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

summary(LM_mean_diff) # get MSE value

eval_models <- as.data.frame(cbind(diff_base_temp$ch4_diff, predict(LM_mean_diff)))

eval_LM_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_mean <- rbind(LM_mean_diff$coefficients)

coefficients_low <- cbind(LM_mean_diff$coefficients[1]-summary(LM_mean_diff)$coefficients[,2][1],
                          LM_mean_diff$coefficients[2]-summary(LM_mean_diff)$coefficients[,2][2])

coefficients_high <- cbind(LM_mean_diff$coefficients[1]+summary(LM_mean_diff)$coefficients[,2][1],
                           LM_mean_diff$coefficients[2]+summary(LM_mean_diff)$coefficients[,2][2])


LM_coefficient_mean <- function(x, p){
  predicted_diff_rate = coefficients_mean[1,1] + p*coefficients_mean[1,2]
  return(predicted_diff_rate)
}

LM_coefficient_low <- function(x, p){
  predicted_diff_rate = coefficients_low[1,1] + p*coefficients_low[1,2]
  return(predicted_diff_rate)
}

LM_coefficient_high <- function(x, p){
  predicted_diff_rate = coefficients_high[1,1] + p*coefficients_high[1,2]
  return(predicted_diff_rate)
}







diff_base_temp <- error_diff %>%
  ungroup(.) %>%
  select(ch4_diff, sd_time, waterbody_type, sd_space, temp_for_model_C) %>%
  mutate(diff_time_low = ifelse(ch4_diff-sd_time < 0, 0.01, ch4_diff-sd_time)) %>%
  mutate(diff_time_high = ifelse(ch4_diff+sd_time < 0, 0.01, ch4_diff+sd_time)) %>%
  mutate(diff_space_low = ifelse(ch4_diff-sd_space < 0, 0.01, ch4_diff-sd_space)) %>%
  mutate(diff_space_high = ifelse(ch4_diff+sd_space < 0, 0.01, ch4_diff+sd_space)) %>%
  na.omit(.)%>%
  filter(waterbody_type == "reservoir")

# First-order Arhennius
mean_diff_FOA_res = nlsLM(ch4_diff ~ A * exp(a * temp_for_model_C),
                      start = list(A = 0.023, a = 0.124),
                      data = diff_base_temp,
                      control = nls.lm.control(maxiter=1000))


summary(mean_diff_FOA_res) # get MSE value

eval_models <- as.data.frame(cbind(diff_base_temp$ch4_diff, predict(mean_diff_FOA_res)))

eval_FOA_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_main <- rbind(mean_diff_FOA_res$m$getPars())


coefficients_low <- cbind(mean_diff_FOA_res$m$getPars()[1]-summary(mean_diff_FOA_res)$coefficients[,2][1],
                          mean_diff_FOA_res$m$getPars()[2]-summary(mean_diff_FOA_res)$coefficients[,2][2])

coefficients_high <- cbind(mean_diff_FOA_res$m$getPars()[1]+summary(mean_diff_FOA_res)$coefficients[,2][1],
                           mean_diff_FOA_res$m$getPars()[2]+summary(mean_diff_FOA_res)$coefficients[,2][2])

FOA_coefficient_mean_res <- function(x){
  predicted_diff_rate = (coefficients_main[1,1] * exp(coefficients_main[1,2]*x))
  return(predicted_diff_rate)
}

FOA_coefficient_low_res <- function(x){
  predicted_diff_rate = (coefficients_main[1,1] * exp(coefficients_main[1,2]*x))
  return(predicted_diff_rate)
}

FOA_coefficient_high_res <- function(x){
  predicted_diff_rate = (coefficients_main[1,1] * exp(coefficients_main[1,2]*x))
  return(predicted_diff_rate)
}


diff_base_temp <- error_diff %>%
  ungroup(.) %>%
  select(ch4_diff, sd_time, waterbody_type, sd_space, surf_area_k) %>%
  mutate(diff_time_low = ifelse(ch4_diff-sd_time < 0, 0.01, ch4_diff-sd_time)) %>%
  mutate(diff_time_high = ifelse(ch4_diff+sd_time < 0, 0.01, ch4_diff+sd_time)) %>%
  mutate(diff_space_low = ifelse(ch4_diff-sd_space < 0, 0.01, ch4_diff-sd_space)) %>%
  mutate(diff_space_high = ifelse(ch4_diff+sd_space < 0, 0.01, ch4_diff+sd_space)) %>%
  na.omit(.)%>%
  filter(waterbody_type == "reservoir")


LM_mean_diff_res= glm(ch4_diff ~ surf_area_k, data = diff_base_temp)

summary(LM_mean_diff_res) # get MSE value

eval_models <- as.data.frame(cbind(diff_base_temp$ch4_diff, predict(LM_mean_diff_res)))

eval_LM_models <- eval_models %>%
  summarize(NSE_mean = NSE(V2, V1))

coefficients_mean <- rbind(LM_mean_diff_res$coefficients)

coefficients_low <- cbind(LM_mean_diff_res$coefficients[1]-summary(LM_mean_diff_res)$coefficients[,2][1],
                          LM_mean_diff_res$coefficients[2]-summary(LM_mean_diff_res)$coefficients[,2][2])

coefficients_high <- cbind(LM_mean_diff_res$coefficients[1]+summary(LM_mean_diff_res)$coefficients[,2][1],
                           LM_mean_diff_res$coefficients[2]+summary(LM_mean_diff_res)$coefficients[,2][2])


LM_coefficient_mean_res <- function(x, p){
  predicted_diff_rate = coefficients_mean[1,1] + p*coefficients_mean[1,2]
  return(predicted_diff_rate)
}

LM_coefficient_low_res <- function(x, p){
  predicted_diff_rate = coefficients_low[1,1] + p*coefficients_low[1,2]
  return(predicted_diff_rate)
}

LM_coefficient_high_res <- function(x, p){
  predicted_diff_rate = coefficients_high[1,1] + p*coefficients_high[1,2]
  return(predicted_diff_rate)
}