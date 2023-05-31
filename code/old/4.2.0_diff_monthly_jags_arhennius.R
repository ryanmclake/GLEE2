ahrennius_diff_model = ("arhennius_diff_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   A ~ dnorm(0,1e-6)
   a ~ dnorm(0,1e-6)
   sd.pro ~ dunif(0, 1000)
   
   #end priors===============================================
   
   for(i in 1:N) {
     
      #process model=============================================
     
      tau.pro[i] <- 1/((sd.pro)*(sd.pro))
      predX[i] <- A*exp(a*(temp[i]))
      X[i] ~ dnorm(predX[i],tau.pro[i])
     
      #end of process model======================================
     
      #data model================================================
     
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
            #end of data model=========================================
   }
  }", file = ahrennius_diff_model)


################################

### MONTHLY TIMESTEP ALL WB ###
diff_base_temp <- base %>% select(lat, lon, ch4_diff, temp_for_model_K, waterbody_id, year, month, tot_sampling_events) %>%
  na.omit(.) %>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  mutate(diff_sd = ifelse(tot_sampling_events < 5,93.6, NA),
         diff_sd = ifelse(tot_sampling_events >= 5,25.3, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 10,2.83, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 15,6.78, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 20,0.819, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 25,0.902, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 30,0.973, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 40,0.0661, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 50,0.00665, diff_sd))


jags.data = list(Y = diff_base_temp$ch4_diff,
                 tau.obs = 1/((diff_base_temp$diff_sd)) ^ 2,
                 N = nrow(diff_base_temp),
                 temp = as.numeric(diff_base_temp$temp_for_model_C))

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = runif(1, 1, 3),
                    A = runif(1, 0.24,1),
                    a = runif(1, 0.023,1),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = ahrennius_diff_model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_diff  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro","A", "a"),
                          n.iter = 20000, n.burnin = 2000, thin = 20)

plot(eval_diff)


all_coefficients <- eval_diff %>%
  spread_draws(sd.pro, A, a) %>%
  filter(.chain == 1)

all_coefficients_DF <- all_coefficients %>%
  summarise(mean_E20 = mean(A),
            sd_E20 = sd(A),
            upper_95_E20 = quantile(A, 0.95, na.rm = T),
            lower_95_E20 = quantile(A, 0.05, na.rm = T),
            mean_omega = mean(a),
            sd_omega = sd(a),
            upper_95_omega = quantile(a, 0.95, na.rm = T),
            lower_95_omega = quantile(a, 0.05, na.rm = T),
            mean_process = mean(sd.pro),
            sd_process = sd(sd.pro),
            upper_95_process = quantile(sd.pro, 0.95, na.rm = T),
            lower_95_process = quantile(sd.pro, 0.05, na.rm = T))


### MONTHLY TIMESTEP JUST LAKES ###
ebu_diff_temp <- base %>% select(waterbody_type, lat, lon, ch4_diff, temp_for_model_K, waterbody_id, year, month, tot_sampling_events) %>%
  na.omit(.) %>%
  filter(waterbody_type == "lake") %>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  mutate(diff_sd = ifelse(tot_sampling_events < 5,93.6, NA),
         diff_sd = ifelse(tot_sampling_events >= 5,25.3, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 10,2.83, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 15,6.78, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 20,0.819, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 25,0.902, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 30,0.973, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 40,0.0661, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 50,0.00665, diff_sd))


jags.data = list(Y = ebu_diff_temp$ch4_diff,
                 tau.obs = 1/((ebu_diff_temp$diff_sd)) ^ 2,
                 N = nrow(ebu_diff_temp),
                 temp = as.numeric(ebu_diff_temp$temp_for_model_C))

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = runif(1, 1, 3),
                    A = runif(1, 0.24,1),
                    a = runif(1, 0.023,1),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = ahrennius_diff_model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_diff  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro","A","a"),
                          n.iter = 20000, n.burnin = 2000, thin = 20)
plot(eval_diff)

lake_coefficients <- eval_diff %>%
  spread_draws(sd.pro, A, a) %>%
  filter(.chain == 1)

lake_coefficients_DF <- lake_coefficients %>%
  summarise(mean_E20 = mean(A),
            sd_E20 = sd(A),
            upper_95_E20 = quantile(A, 0.95, na.rm = T),
            lower_95_E20 = quantile(A, 0.05, na.rm = T),
            mean_omega = mean(a),
            sd_omega = sd(a),
            upper_95_omega = quantile(a, 0.95, na.rm = T),
            lower_95_omega = quantile(a, 0.05, na.rm = T),
            mean_process = mean(sd.pro),
            sd_process = sd(sd.pro),
            upper_95_process = quantile(sd.pro, 0.95, na.rm = T),
            lower_95_process = quantile(sd.pro, 0.05, na.rm = T))



### MONTHLY TIMESTEP JUST RESERVOIRS ###
ebu_diff_temp <- base %>% select(waterbody_type, lat, lon, ch4_diff, temp_for_model_K, waterbody_id, year, month, tot_sampling_events) %>%
  na.omit(.) %>%
  filter(waterbody_type == "reservoir") %>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  mutate(diff_sd = ifelse(tot_sampling_events < 5,93.6, NA),
         diff_sd = ifelse(tot_sampling_events >= 5,25.3, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 10,2.83, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 15,6.78, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 20,0.819, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 25,0.902, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 30,0.973, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 40,0.0661, diff_sd),
         diff_sd = ifelse(tot_sampling_events > 50,0.00665, diff_sd))


jags.data = list(Y = ebu_diff_temp$ch4_diff,
                 tau.obs = 1/((ebu_diff_temp$diff_sd)) ^ 2,
                 N = nrow(ebu_diff_temp),
                 temp = as.numeric(ebu_diff_temp$temp_for_model_C))

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = runif(1, 1, 3),
                    A = runif(1, 0.24,1),
                    a = runif(1, 0.023,1),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = ahrennius_diff_model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_diff  <- coda.samples(model = j.model,
                           variable.names = c("sd.pro","A","a"),
                           n.iter = 200000, n.burnin = 20000, thin = 200)
plot(eval_diff)

res_coefficients <- eval_diff %>%
  spread_draws(sd.pro, A, a) %>%
  filter(.chain == 1)

res_coefficients_DF <- res_coefficients %>%
  summarise(mean_E20 = mean(A),
            sd_E20 = sd(A),
            upper_95_E20 = quantile(A, 0.95, na.rm = T),
            lower_95_E20 = quantile(A, 0.05, na.rm = T),
            mean_omega = mean(a),
            sd_omega = sd(a),
            upper_95_omega = quantile(a, 0.95, na.rm = T),
            lower_95_omega = quantile(a, 0.05, na.rm = T),
            mean_process = mean(sd.pro),
            sd_process = sd(sd.pro),
            upper_95_process = quantile(sd.pro, 0.95, na.rm = T),
            lower_95_process = quantile(sd.pro, 0.05, na.rm = T))
