ahrennius_ebu_model = ("arhennius_ebu_model.txt")
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
      predX[i] <- A * (a^(temp[i]-20))
      X[i] ~ dnorm(predX[i],tau.pro[i])
     
      #end of process model======================================
     
      #data model================================================
     
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
            #end of data model=========================================
   }
  }", file = ahrennius_ebu_model)

### MONTHLY TIMESTEP ###
ebu_base_temp <- base %>% select(ch4_ebu, temp_for_model_K, waterbody_id, month, tot_sampling_events) %>%
  na.omit(.) %>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events <= 10, 29.434313, NA)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 10, 8.903169, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 20, 1.914520, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 30, 1.131719, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 40, 1.123484, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 50, 1.101714, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 60, 1.154789, ebu_sd)) %>%
  mutate(ebu_sd = ifelse(tot_sampling_events > 100, 1.151526, ebu_sd)) %>%
  mutate(sd2 = 0.001)


jags.data = list(Y = ebu_base_temp$ch4_ebu,
                 tau.obs = 1/((ebu_base_temp$sd2)) ^ 2,
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
print(gelman.diag(eval_ebu))


parameter <- eval_ebu %>%
  spread_draws(sd.pro, A, a) %>%
  select(sd.pro, A, a)



ebu_arhennius_jags <- function(A, a, temp, Q){
  est = A * (a^(temp-20)) + rnorm(1000,0, sd = Q)
  return(est)
}

parms <- sample_n(parameter, 1000, replace=TRUE)

temps <- c(seq(from = -10, to = 35, by = 1))

out <- list()

for(s in 1:length(temps)){

  prediction <- var(ebu_arhennius_jags(temp = temps[s],
                                       A = mean(parms$A),
                                       a = mean(parms$a),
                                       Q = parms$sd.pro))
  out[[s]] <- validation
}



model_validate = as.data.frame(do.call(rbind, out)) %>%
  pivot_longer(!`ebu_base_temp$temp_for_model_C`, names_to = "iteration", values_to = "value") %>%
  rename(temp = `ebu_base_temp$temp_for_model_C`) %>%
  ungroup() %>%
  select(-iteration) %>%
  group_by(temp) %>%
  summarize(mean_eb_flux = mean(value),
            sd_eb_flux = sd(value), 
            var_flux = var(value))




### MONTHLY TIMESTEP ###
ebu_base_temp <- base %>% select(ch4_ebu, temp_for_model_K, waterbody_id, month, tot_sampling_events) %>%
  na.omit(.) %>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  mutate(ebu_sd = 0.01)


jags.data = list(Y = ebu_base_temp$ch4_ebu,
                 tau.obs = 1/((ebu_base_temp$ebu_sd)) ^ 2,
                 N = nrow(ebu_base_temp),
                 temp = ebu_base_temp$temp_for_model_C)

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = runif(1, 0.01, 2),
                    A = runif(1, 4.14072,10),
                    a = runif(1, 4.895e-02, 3),
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
print(gelman.diag(eval_ebu))


parameter <- eval_ebu %>%
  spread_draws(sd.pro, A, a) %>%
  select(sd.pro, A, a)

parms <- sample_n(parameter, 1000, replace=TRUE)

ebu_arhennius_jags <- function(A, a, temp, Q){
  est = (A * exp(a*temp)) + rnorm(1000,0, sd = Q)
  return(est)
}

out <- list()

for(s in 1:length(ebu_base_temp$temp_for_model_C)){
  
  validation <- ebu_arhennius_jags(temp = ebu_base_temp$temp_for_model_C,
                                   A = parms$A,
                                   a = parms$a,
                                   Q = parms$sd.pro)
  out[[s]] <- validation
}

model_validate = as.data.frame(do.call(rbind, out)) %>%
  cbind(., ebu_base_temp$temp_for_model_C) %>%
  pivot_longer(!`ebu_base_temp$temp_for_model_C`, names_to = "iteration", values_to = "value") %>%
  rename(temp = `ebu_base_temp$temp_for_model_C`) %>%
  ungroup() %>%
  select(-iteration) %>%
  group_by(temp) %>%
  summarize(mean_eb_flux = mean(value),
            sd_eb_flux = sd(value), 
            var_flux = var(value))
