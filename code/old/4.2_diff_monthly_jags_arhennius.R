ahrennius_ebu_model = ("arhennius_ebu_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   A ~ dnorm(0,1e-6)
   a ~ dnorm(0,1e-6)
   #sd.pro ~ dgamma(0.01, 0.01)
   sd.pro ~ dunif(0, 1000)
   
   #end priors===============================================
   
   for(i in 1:N) {
     
      #process model=============================================
     
      tau.pro[i] <- 1/((sd.pro)*(sd.pro))
      #log(predX[i]) <- A+a*(temp[i])
      predX[i] <- A*exp(a*(temp[i]))
      X[i] ~ dnorm(predX[i],tau.pro[i])
     
      #end of process model======================================
     
      #data model================================================
     
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
            #end of data model=========================================
   }
  }", file = ahrennius_ebu_model)


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

j.model   <- jags.model(file = ahrennius_ebu_model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_diff  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro","A", "a"),
                          n.iter = 200000, n.burnin = 20000, thin = 200)



plot(eval_diff)
print(gelman.diag(eval_diff))

parameter <- eval_diff %>%
  spread_draws(sd.pro, A, a) %>%
  select(sd.pro, A, a)

summary(parameter)

sd(parameter$A)

diff_arhennius <- function(temp){
  est = (rnorm(1,mean(parameter$A), sd = sd(parameter$A)) * exp(rnorm(1,mean(parameter$a), sd = sd(parameter$a))*temp)) + rnorm(1,0, sd = sd(parameter$sd.pro))
  return(est)
}

out_jags <- list()
iteration <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)

for(i in 1:length(iteration)){
  
  d <- diff_base_temp %>%
    mutate(ebu_model_prediction = mapply(diff_arhennius, temp_for_model_C),
           'date' = lubridate::make_date(year = year, month = month),
           iteration = iteration[i])
  
  out_jags[[i]] <- d
  
}

dat <- as.data.frame(do.call(rbind, out_jags)) %>%
  group_by(lat, lon) %>%
  summarize(mean_diff_prediction = mean(ebu_model_prediction),
            sd_diff_prediction = sd(ebu_model_prediction))

validate <- left_join(diff_base_temp, dat, by = c("lat", "lon")) %>%
  summarise(NSE = NSE(mean_diff_prediction, ch4_diff),
            rmse = rmse(mean_diff_prediction, ch4_diff))


### MONTHLY TIMESTEP JUST LAKES ###
ebu_base_temp <- base %>% select(waterbody_type, lat, lon, ch4_ebu, temp_for_model_K, waterbody_id, year, month, tot_sampling_events) %>%
  na.omit(.) %>%
  filter(waterbody_type == "lake") %>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events < 5,388, NA),
         ebu_sd = ifelse(tot_sampling_events >= 5,43.1, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 10,35.3, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 15,5.61, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 20,4.26, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 25,0.902, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 30,0.248, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 50,0.0426, ebu_sd)) 


jags.data = list(Y = ebu_base_temp$ch4_ebu,
                 tau.obs = 1/((ebu_base_temp$ebu_sd)) ^ 2,
                 N = nrow(ebu_base_temp),
                 temp = as.numeric(ebu_base_temp$temp_for_model_C))

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = runif(1, 0.01, 2),
                    A = runif(1, 100,150),
                    a = runif(1, 1, 3),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = ahrennius_ebu_model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_ebu  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro","A","a"),
                          n.iter = 500000, n.burnin = 50000, thin = 500)
plot(eval_ebu)

parameter <- eval_ebu %>%
  spread_draws(sd.pro, A, a) %>%
  select(sd.pro, A, a)

ebu_arhennius <- function(temp){
  est = rnorm(1,mean(parameter$A), sd = sd(parameter$A)) * rnorm(1,mean(parameter$a), sd = sd(parameter$a)) ^ (temp-20) + rnorm(1,0, sd = 8.752777)
  return(est)
}

out_jags <- list()
iteration <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)

for(i in 1:length(iteration)){
  
  d <- ebu_base_temp %>%
    mutate(ebu_model_prediction = mapply(ebu_arhennius, temp_for_model_C),
           'date' = lubridate::make_date(year = year, month = month),
           iteration = iteration[i])
  
  out_jags[[i]] <- d
  
}

dat <- as.data.frame(do.call(rbind, out_jags)) %>%
  group_by(lat, lon) %>%
  summarize(mean_ebu_prediction = mean(ebu_model_prediction),
            sd_ebu_prediction = sd(ebu_model_prediction))

validate <- left_join(ebu_base_temp, dat, by = c("lat", "lon")) %>%
  summarise(NSE = NSE(mean_ebu_prediction, ch4_ebu),
            rmse = rmse(mean_ebu_prediction, ch4_ebu))



### MONTHLY TIMESTEP JUST RESERVOIRS ###
ebu_base_temp <- base %>% select(waterbody_type, lat, lon, ch4_ebu, temp_for_model_K, waterbody_id, year, month, tot_sampling_events) %>%
  na.omit(.) %>%
  filter(waterbody_type == "reservoir") %>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events < 5,388, NA),
         ebu_sd = ifelse(tot_sampling_events >= 5,43.1, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 10,35.3, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 15,5.61, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 20,4.26, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 25,0.902, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 30,0.248, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 50,0.0426, ebu_sd)) 


jags.data = list(Y = ebu_base_temp$ch4_ebu,
                 tau.obs = 1/((ebu_base_temp$ebu_sd)) ^ 2,
                 N = nrow(ebu_base_temp),
                 temp = as.numeric(ebu_base_temp$temp_for_model_C))

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = runif(1, 0.01, 2),
                    A = runif(1, 100,150),
                    a = runif(1, 1, 3),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = ahrennius_ebu_model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_ebu  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro","A","a"),
                          n.iter = 200000, n.burnin = 20000, thin = 200)
plot(eval_ebu)

parameter <- eval_ebu %>%
  spread_draws(sd.pro, A, a) %>%
  select(sd.pro, A, a)

ebu_arhennius <- function(temp){
  est = rnorm(1,mean(parameter$A), sd = sd(parameter$A)) * rnorm(1,mean(parameter$a), sd = sd(parameter$a)) ^ (temp-20) + rnorm(1,0, sd = 8.752777)
  return(est)
}

out_jags <- list()
iteration <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)

for(i in 1:length(iteration)){
  
  d <- ebu_base_temp %>%
    mutate(ebu_model_prediction = mapply(ebu_arhennius, temp_for_model_C),
           'date' = lubridate::make_date(year = year, month = month),
           iteration = iteration[i])
  
  out_jags[[i]] <- d
  
}

dat <- as.data.frame(do.call(rbind, out_jags)) %>%
  group_by(lat, lon) %>%
  summarize(mean_ebu_prediction = mean(ebu_model_prediction),
            sd_ebu_prediction = sd(ebu_model_prediction))

validate <- left_join(ebu_base_temp, dat, by = c("lat", "lon")) %>%
  summarise(NSE = NSE(mean_ebu_prediction, ch4_ebu),
            rmse = rmse(mean_ebu_prediction, ch4_ebu))

################################

################################

### YEARLY TIMESTEP ALL WB ###
ebu_base_temp <- base %>% select(lat, lon, ch4_ebu, temp_for_model_K, waterbody_id, year, month, tot_sampling_events) %>%
  na.omit(.) %>%
  group_by(lat, lon, year) %>%
  summarize(tot_sampling_events = sum(tot_sampling_events),
            temp_for_model_K = mean(temp_for_model_K),
            ch4_ebu = mean(ch4_ebu)) %>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events < 5,388, NA),
         ebu_sd = ifelse(tot_sampling_events >= 5,43.1, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 10,35.3, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 15,5.61, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 20,4.26, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 25,0.902, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 30,0.248, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 50,0.0426, ebu_sd)) 


jags.data = list(Y = ebu_base_temp$ch4_ebu,
                 tau.obs = 1/((ebu_base_temp$ebu_sd)) ^ 2,
                 N = nrow(ebu_base_temp),
                 temp = as.numeric(ebu_base_temp$temp_for_model_C))

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = runif(1, 0.01, 2),
                    A = runif(1, 100,150),
                    a = runif(1, 1, 3),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = ahrennius_ebu_model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_ebu  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro","A","a"),
                          n.iter = 200000, n.burnin = 20000, thin = 200)
plot(eval_ebu)

parameter <- eval_ebu %>%
  spread_draws(sd.pro, A, a) %>%
  select(sd.pro, A, a)

ebu_arhennius <- function(temp){
  est = rnorm(1,mean(parameter$A), sd = sd(parameter$A)) * rnorm(1,mean(parameter$a), sd = sd(parameter$a)) ^ (temp-20) + rnorm(1,0, sd = 8.752777)
  return(est)
}

out_jags <- list()
iteration <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)

for(i in 1:length(iteration)){
  
  d <- ebu_base_temp %>%
    mutate(ebu_model_prediction = mapply(ebu_arhennius, temp_for_model_C),
           year = year,
           iteration = iteration[i])
  
  out_jags[[i]] <- d
  
}

dat <- as.data.frame(do.call(rbind, out_jags)) %>%
  group_by(lat, lon) %>%
  summarize(mean_ebu_prediction = mean(ebu_model_prediction, na.rm = T),
            sd_ebu_prediction = sd(ebu_model_prediction))

validate <- left_join(ebu_base_temp, dat, by = c("lat", "lon")) %>%
  ungroup(.) %>%
  summarise(NSE = NSE(mean_ebu_prediction, ch4_ebu),
            rmse = rmse(mean_ebu_prediction, ch4_ebu))


### YEARLY TIMESTEP LAKES ###
ebu_base_temp <- base %>% select(lat, lon, ch4_ebu, waterbody_type, temp_for_model_K, waterbody_id, year, month, tot_sampling_events) %>%
  na.omit(.) %>%
  filter(waterbody_type == "lake") %>%
  group_by(lat, lon, year) %>%
  summarize(tot_sampling_events = sum(tot_sampling_events),
            temp_for_model_K = mean(temp_for_model_K),
            ch4_ebu = mean(ch4_ebu)) %>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events < 5,388, NA),
         ebu_sd = ifelse(tot_sampling_events >= 5,43.1, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 10,35.3, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 15,5.61, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 20,4.26, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 25,0.902, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 30,0.248, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 50,0.0426, ebu_sd)) 


jags.data = list(Y = ebu_base_temp$ch4_ebu,
                 tau.obs = 1/((ebu_base_temp$ebu_sd)) ^ 2,
                 N = nrow(ebu_base_temp),
                 temp = as.numeric(ebu_base_temp$temp_for_model_C))

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = runif(1, 0.01, 2),
                    A = runif(1, 100,150),
                    a = runif(1, 1, 3),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = ahrennius_ebu_model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_ebu  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro","A","a"),
                          n.iter = 200000, n.burnin = 20000, thin = 200)
plot(eval_ebu)

parameter <- eval_ebu %>%
  spread_draws(sd.pro, A, a) %>%
  select(sd.pro, A, a)

ebu_arhennius <- function(temp){
  est = rnorm(1,mean(parameter$A), sd = sd(parameter$A)) * rnorm(1,mean(parameter$a), sd = sd(parameter$a)) ^ (temp-20) + rnorm(1,0, sd = 8.752777)
  return(est)
}

out_jags <- list()
iteration <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)

for(i in 1:length(iteration)){
  
  d <- ebu_base_temp %>%
    mutate(ebu_model_prediction = mapply(ebu_arhennius, temp_for_model_C),
           year = year,
           iteration = iteration[i])
  
  out_jags[[i]] <- d
  
}

dat <- as.data.frame(do.call(rbind, out_jags)) %>%
  group_by(lat, lon) %>%
  summarize(mean_ebu_prediction = mean(ebu_model_prediction, na.rm = T),
            sd_ebu_prediction = sd(ebu_model_prediction))

validate <- left_join(ebu_base_temp, dat, by = c("lat", "lon")) %>%
  ungroup(.) %>%
  summarise(NSE = NSE(mean_ebu_prediction, ch4_ebu),
            rmse = rmse(mean_ebu_prediction, ch4_ebu))



### MONTHLY TIMESTEP JUST RESERVOIRS ###
ebu_base_temp <- base %>% select(waterbody_type, lat, lon, ch4_ebu, temp_for_model_K, waterbody_id, year, month, tot_sampling_events) %>%
  na.omit(.) %>%
  filter(waterbody_type == "reservoir") %>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events < 5,388, NA),
         ebu_sd = ifelse(tot_sampling_events >= 5,43.1, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 10,35.3, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 15,5.61, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 20,4.26, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 25,0.902, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 30,0.248, ebu_sd),
         ebu_sd = ifelse(tot_sampling_events > 50,0.0426, ebu_sd)) 


jags.data = list(Y = ebu_base_temp$ch4_ebu,
                 tau.obs = 1/((ebu_base_temp$ebu_sd)) ^ 2,
                 N = nrow(ebu_base_temp),
                 temp = as.numeric(ebu_base_temp$temp_for_model_C))

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = runif(1, 0.01, 2),
                    A = runif(1, 100,150),
                    a = runif(1, 1, 3),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = ahrennius_ebu_model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_ebu  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro","A","a"),
                          n.iter = 200000, n.burnin = 20000, thin = 200)
plot(eval_ebu)

parameter <- eval_ebu %>%
  spread_draws(sd.pro, A, a) %>%
  select(sd.pro, A, a)

ebu_arhennius <- function(temp){
  est = rnorm(1,mean(parameter$A), sd = sd(parameter$A)) * rnorm(1,mean(parameter$a), sd = sd(parameter$a)) ^ (temp-20) + rnorm(1,0, sd = 8.752777)
  return(est)
}

out_jags <- list()
iteration <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)

for(i in 1:length(iteration)){
  
  d <- ebu_base_temp %>%
    mutate(ebu_model_prediction = mapply(ebu_arhennius, temp_for_model_C),
           'date' = lubridate::make_date(year = year, month = month),
           iteration = iteration[i])
  
  out_jags[[i]] <- d
  
}

dat <- as.data.frame(do.call(rbind, out_jags)) %>%
  group_by(lat, lon) %>%
  summarize(mean_ebu_prediction = mean(ebu_model_prediction),
            sd_ebu_prediction = sd(ebu_model_prediction))

validate <- left_join(ebu_base_temp, dat, by = c("lat", "lon")) %>%
  summarise(NSE = NSE(mean_ebu_prediction, ch4_ebu),
            rmse = rmse(mean_ebu_prediction, ch4_ebu))