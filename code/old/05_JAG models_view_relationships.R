
if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, MCMCvis, lubridate, tidybayes,
               ncdf4, reshape2, zoo, patchwork, hydroGOF, viridis,
               imputeTS, devtools, scales, forecast, coda, rjags, R2jags,gridExtra)

w <- read.table("./data/organized_data_to_append/global_lake_res_DB_refs_WWF_GLCP.txt", header = TRUE, row.names = NULL) %>%
  select(ref, mean_diff, mean_ebu, ebu_sd, diff_sd, mean_annual_temp_k, mean_obs_wtemp_k, elevation, lat, waterbody_type, country, biome_type, hybas_id) %>%
  mutate(temp_for_model = ifelse(is.na(mean_annual_temp_k),mean_obs_wtemp_k, mean_annual_temp_k)) 


ebu <- w %>% select(country, biome_type, mean_ebu, ebu_sd, temp_for_model, waterbody_type, lat) %>% na.omit(.) %>% filter(mean_ebu >= 0) %>%
  mutate(temp_for_model = temp_for_model - 273.15)%>%
  mutate(temp_for_model = ifelse(temp_for_model < -100, NA, temp_for_model)) %>%
  na.omit(.)

diff <- w %>% select(country, biome_type, mean_diff, diff_sd, temp_for_model, lat) %>% na.omit(.) %>% filter(mean_diff >= 0) %>%
  mutate(temp_for_model = temp_for_model - 273.15)%>%
  mutate(temp_for_model = ifelse(temp_for_model < -100, NA, temp_for_model)) %>%
  na.omit(.) %>%
  filter(mean_diff < 5000)



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
                 tau.obs = 1/((ebu$ebu_sd+.1)) ^ 2,
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
temp_ebu_function_mcclure_mean <- function(x) {(38.39 * 1.0133 ^ (x - 20))}
temp_ebu_function_mcclure_min <- function(x) {(36.77 * 1.0100 ^ (x - 20))}
temp_ebu_function_mcclure_max <- function(x) {(39.91 * 1.0166 ^ (x - 20))}

ebu_function_output <- data.frame(x = -20:40,            # Create data for ggplot2
                                  values = c(temp_ebu_function_johnson(-20:40),
                                             temp_ebu_function_mcclure_mean(-20:40),
                                             temp_ebu_function_mcclure_min(-20:40),
                                             temp_ebu_function_mcclure_max(-20:40)),
                                  model = rep(c("Johnson - No Bayes",
                                                "Bayes - Mean Parameter",
                                                "Bayes - 1st Q Parameters",
                                                "Bayes - 3rd Q Parameter"), each = 61))

ggplot(ebu) +
  geom_point(aes(temp_for_model, mean_ebu)) +
  geom_line(data = ebu_function_output, aes(x, values, group = model, color = model), lwd = 2)+
  scale_color_viridis(discrete = T, option = "C")+
  theme_classic()




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

parameter_diff <- eval_diff %>%
  spread_draws(sd.pro, alpha, beta) %>%
  select(sd.pro, alpha, beta)

summary(parameter_diff)


temp_diff_function_johnson <- function(x) {0.023*exp(0.124*x)}
temp_diff_function_mcclure_mean <- function(x) {2.66978*exp(0.11747*x)}
temp_diff_function_mcclure_min <- function(x) {0.48161*exp(0.07117*x)}
temp_diff_function_mcclure_max <- function(x) {4.78967*exp(0.16094*x)}

diff_function_output <- data.frame(x = -20:40,            # Create data for ggplot2
                                   values = c(temp_diff_function_johnson(-20:40),
                                              temp_diff_function_mcclure_mean(-20:40),
                                              temp_diff_function_mcclure_min(-20:40),
                                              temp_diff_function_mcclure_max(-20:40)),
                                   model = rep(c("Matt Johnson - No Bayes",
                                                 "Bayes - Mean Parameter",
                                                 "Bayes - 1st Q Parameter",
                                                 "Bayes - 3rd Q Parameter"), each = 61))
ggplot(diff) +
  geom_point(aes(temp_for_model, mean_diff)) +
  geom_line(data = diff_function_output, aes(x, values, group = model, color = model), lwd = 2)+
  scale_color_viridis(discrete = T, option = "D")+
  theme_classic()+
  ylim(c(0,1500))


