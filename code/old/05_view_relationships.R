
if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, MCMCvis, lubridate, tidybayes,
               ncdf4, reshape2, zoo, patchwork, hydroGOF, viridis,
               imputeTS, devtools, scales, forecast, coda, rjags, R2jags,gridExtra)

w <- read.table("./data/organized_data_to_append/global_lake_res_DB_refs_WWF_GLCP.txt", header = TRUE, row.names = NULL) %>%
  select(mean_diff, mean_ebu, ebu_sd, diff_sd, mean_annual_temp_k, mean_obs_wtemp_k, lat, country, biome_type) %>%
  mutate(temp_for_model = ifelse(is.na(mean_annual_temp_k),mean_obs_wtemp_k, mean_annual_temp_k)) 


ebu <- w %>% select(country, biome_type, mean_ebu, ebu_sd, temp_for_model, lat) %>% na.omit(.) %>% filter(mean_ebu >= 0) %>%
  mutate(temp_for_model = temp_for_model - 273.15)%>%
  mutate(temp_for_model = ifelse(temp_for_model < -100, NA, temp_for_model)) %>%
  na.omit(.)

diff <- w %>% select(country, biome_type, mean_diff, diff_sd, temp_for_model, lat) %>% na.omit(.) %>% filter(mean_diff >= 0) %>%
  mutate(temp_for_model = temp_for_model - 273.15)%>%
  mutate(temp_for_model = ifelse(temp_for_model < -100, NA, temp_for_model)) %>%
  na.omit(.)



model= ("ebu_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   ebtwenty ~ dunif(0,1000)
   theta ~ dunif(0,100)
   sd.pro ~ dunif(0, 1000)
   
   #end priors===============================================
   
   for(i in 1:N) {
     
      #process model=============================================
     
      tau.pro[i] <- 1/((sd.pro)*(sd.pro))
      predX[i] <- ebtwenty * (theta^(temp[i]-20))
      X[i] ~ dnorm(predX[i],tau.pro[i])
     
      #end of process model======================================
     
      #data model================================================
     
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
            #end of data model=========================================
   }
  }", file = model)


set.seed(1)

# cal <- eb_mid1 %>% sample_n(., 100)
# cal_vals <- c(cal$row_names)
# val <- eb_mid1 %>% filter(!row_names %in% cal_vals) %>%
#   select(-row_names) %>%
#   tibble::rownames_to_column(., "row_names") %>%
#   mutate(row_names = as.numeric(row_names))


jags.data = list(Y = ebu$mean_ebu,
                 tau.obs = 1/((ebu$ebu_sd+1)) ^ 2,
                 N = nrow(ebu),
                 temp = ebu$temp_for_model)

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = runif(1, 0.01, 2),
                    theta = runif(1, 1.2,1.4),
                    ebtwenty = runif(1, 50, 150),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_ebu  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro","ebtwenty","theta"),
                          n.iter = 200000, n.burnin = 20000, thin = 200)

plot(eval_ebu)
print(gelman.diag(eval_ebu))

parameter <- eval_ebu %>%
  spread_draws(sd.pro, ebtwenty, theta) %>%
  select(sd.pro, ebtwenty, theta)

summary(parameter)



temp_ebu_function_johnson <- function(x) {(100 * 1.1 ^ (x - 20))}
temp_ebu_function_mcclure_mean <- function(x) {(39.42 * 1.0107 ^ (x - 20))}
temp_ebu_function_mcclure_min <- function(x) {(28.33 * 0.9948 ^ (x - 20))}
temp_ebu_function_mcclure_max <- function(x) {(51.48 * 1.0283 ^ (x - 20))}

ebu_function_output <- data.frame(x = -10:40,            # Create data for ggplot2
                                  values = c(temp_ebu_function_johnson(-10:40),
                                             temp_ebu_function_mcclure_mean(-10:40),
                                             temp_ebu_function_mcclure_min(-10:40),
                                             temp_ebu_function_mcclure_max(-10:40)),
                                  model = rep(c("Johnson",
                                                "Bayes - Mean",
                                                "Bayes - Min",
                                                "Bayes - Max"), each = 51))

ggplot(ebu) +
  geom_point(aes(temp_for_model, mean_ebu)) +
  geom_line(data = ebu_function_output, aes(x, values, group = model, color = model))+
  scale_color_viridis(discrete = T, option = "C")


model= ("diff_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   alpha ~ dunif(0,10000)
   beta ~ dunif(0,1000)
   sd.pro ~ dunif(0, 1000)
   
   #end priors===============================================
   
   for(i in 1:N) {
     
      #process model=============================================
     
      tau.pro[i] <- 1/((sd.pro)*(sd.pro))
      predX[i] <- alpha*exp(beta*(temp[i]))
      X[i] ~ dnorm(predX[i],tau.pro[i])
     
      #end of process model======================================
     
      #data model================================================
     
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
            #end of data model=========================================
   }
  }", file = model)


set.seed(1)

# cal <- eb_mid1 %>% sample_n(., 100)
# cal_vals <- c(cal$row_names)
# val <- eb_mid1 %>% filter(!row_names %in% cal_vals) %>%
#   select(-row_names) %>%
#   tibble::rownames_to_column(., "row_names") %>%
#   mutate(row_names = as.numeric(row_names))


jags.data = list(Y = diff$mean_diff,
                 tau.obs = 1/((diff$diff_sd+1)) ^ 2,
                 N = nrow(ebu),
                 temp = diff$temp_for_model)

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = runif(1, 0.01, 2),
                    alpha = runif(1, 0.023,2),
                    beta = runif(1, 0.124, 1),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_diff  <- coda.samples(model = j.model,
                           variable.names = c("sd.pro","alpha","beta"),
                           n.iter = 200000, n.burnin = 20000, thin = 200)

plot(eval_diff)
print(gelman.diag(eval_diff))

parameter <- eval_diff %>%
  spread_draws(sd.pro, alpha, beta) %>%
  select(sd.pro, alpha, beta)

summary(parameter)


temp_diff_function_johnson <- function(x) {0.023*exp(0.124*x)}
temp_diff_function_mcclure_mean <- function(x) {7.103*exp(0.04979*x)}
temp_diff_function_mcclure_min <- function(x) {1.726*exp(0.01708*x)}
temp_diff_function_mcclure_max <- function(x) {14.048*exp(0.10563*x)}

diff_function_output <- data.frame(x = -10:40,            # Create data for ggplot2
                                   values = c(temp_diff_function_johnson(-10:40),
                                              temp_diff_function_mcclure_mean(-10:40),
                                              temp_diff_function_mcclure_min(-10:40),
                                              temp_diff_function_mcclure_max(-10:40)),
                                   model = rep(c("Johnson",
                                                 "Bayes - Mean",
                                                 "Bayes - Min",
                                                 "Bayes - Max"), each = 51))
ggplot(diff) +
  geom_point(aes(temp_for_model, mean_diff)) +
  geom_line(data = diff_function_output, aes(x, values, group = model, color = model))+
  scale_color_viridis(discrete = T, option = "C")





# 
# 
# 
# 
# 
# 
# if (!require("pacman")) install.packages("pacman")
# pacman::p_load(tidyverse,viridis)
# 
# temp_ebu_function_johnson <- function(x) {(100 * 1.1 ^ (x - 20))}
# temp_ebu_function_mcclure_mean <- function(x) {(125 * 1.3 ^ (x - 20))}
# temp_ebu_function_mcclure_min <- function(x) {(125 * 1.2 ^ (x - 20))}
# temp_ebu_function_mcclure_max <- function(x) {(125 * 1.4 ^ (x - 20))}
# 
# ebu_function_output <- data.frame(x = -10:40,            # Create data for ggplot2
#                                   values = c(temp_ebu_function_johnson(-10:40),
#                                              temp_ebu_function_mcclure_mean(-10:40),
#                                              temp_ebu_function_mcclure_min(-10:40),
#                                              temp_ebu_function_mcclure_max(-10:40)),
#                                   model = rep(c("Johnson",
#                                                 "Bayes - Mean",
#                                                 "Bayes - Min",
#                                                 "Bayes - Max"), each = 51))
# 
# 
# w <- read.table("./data/organized_data_to_append/global_lake_res_DB_refs_WWF_GLCP.txt", header = TRUE, row.names = NULL) %>%
#   select(mean_diff, mean_ebu, ebu_sd, diff_sd, mean_annual_temp_k, mean_obs_wtemp_k, lat, country, biome_type) %>%
#   mutate(temp_for_model = ifelse(is.na(mean_annual_temp_k),mean_obs_wtemp_k, mean_annual_temp_k)) 
# 
# 
# ebu <- w %>% select(country, biome_type, mean_ebu, ebu_sd, temp_for_model, lat) %>% na.omit(.) %>% filter(mean_ebu >= 0) %>%
#   mutate(temp_for_model = temp_for_model - 273.15)%>%
#   mutate(temp_for_model = ifelse(temp_for_model < -100, NA, temp_for_model)) %>%
#   na.omit(.)
# 
# ggplot(ebu) +
#   geom_point(aes(temp_for_model, mean_ebu)) +
#   geom_line(data = ebu_function_output, aes(x, values, group = model, color = model))+
#   scale_color_viridis(discrete = T, option = "C")+
#   ylim(c(0,7000))
# 
# 
# 
# 
# diff <- w %>% select(mean_diff, diff_sd, temp_for_model, lat) %>% na.omit(.) %>% filter(mean_diff > 0)
# 
# 
# 
# 
# 
# 
# w %>% select(mean_diff, temp_for_model, lat) %>%
#   mutate(temp_for_model = temp_for_model-273.15) %>%
#   filter(mean_diff < 5000) %>%
#   ggplot(.) +
#   geom_point(aes(temp_for_model, mean_diff, color = lat)) +
#   geom_line(data = diff_function_johnson, aes(x, values))
#   
# 
# w %>% select(mean_ebu, temp_for_model, lat) %>%
#   mutate(temp_for_model = temp_for_model-273.15) %>%
#   
# 
# 
# 
# 
# w %>% select(mean_diff, temp_for_model, lat) %>%
#   filter(lat < 45) %>%
#   filter(lat > 0) %>%
#   ggplot(.) +
#   geom_point(aes(temp_for_model, mean_diff, color = lat))
# 
# 
# 
# w %>% select(mean_diff, temp_for_model, lat) %>%
#   filter(lat < 0) %>%
#   ggplot(.) +
#   geom_point(aes(temp_for_model, mean_diff, color = lat))
# 
# 
# 
# 
# 
# 
# w %>% select(mean_diff, temp_for_model, lat) %>%
#   filter(lat < -45) %>%
#   ggplot(.) +
#   geom_point(aes(temp_for_model, mean_diff, color = lat))
