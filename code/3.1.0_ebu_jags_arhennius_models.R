

ebu_coefficients <- read_csv("./data/ebullition_coefficients.csv")
summary(ebu_coefficients)
sd(ebu_coefficients$E20)
sd(ebu_coefficients$omega)

ahrennius_ebu_model = ("arhennius_ebu_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   A ~ dnorm(319.5,1/829.4559)
   a ~ dnorm(1.110,1/0.128368)
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




################################

### MONTHLY TIMESTEP ALL WB ###
ebu_base_temp <- base %>% select(lat, lon, ch4_ebu, temp_for_model_K, waterbody_type, year, month, tot_sampling_events) %>%
  na.omit(.) %>%
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
                          n.iter = 20000, n.burnin = 2000, thin = 20)
plot(eval_ebu)

all_coefficients <- eval_ebu %>%
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










### MONTHLY TIMESTEP LAKES ###
ebu_base_temp <- base %>% select(lat, lon, ch4_ebu, temp_for_model_K, waterbody_type, year, month, tot_sampling_events) %>%
  filter(waterbody_type == "lake") %>%
  na.omit(.) %>%
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
  init[[i]] <- list(sd.pro = runif(1, 1, 2),
                    A = runif(1, mean(ebu_coefficients$E20), sd(ebu_coefficients$E20)),
                    a = runif(1, mean(ebu_coefficients$omega), 3),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = ahrennius_ebu_model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_ebu  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro","A","a"),
                          n.iter = 20000, n.burnin = 2000, thin = 20)
plot(eval_ebu)


lake_coefficients <- eval_ebu %>%
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











### MONTHLY TIMESTEP RESERVOIRS ###
ebu_base_temp <- base %>% select(lat, lon, ch4_ebu, temp_for_model_K, waterbody_type, year, month, tot_sampling_events) %>%
  filter(waterbody_type == "reservoir") %>%
  na.omit(.) %>%
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
  init[[i]] <- list(sd.pro = runif(1, 1, 2),
                    A = runif(1, mean(ebu_coefficients$E20), sd(ebu_coefficients$E20)),
                    a = runif(1, mean(ebu_coefficients$omega), 3),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = ahrennius_ebu_model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_ebu  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro","A","a"),
                          n.iter = 20000, n.burnin = 2000, thin = 20)
plot(eval_ebu)


res_coefficients <- eval_ebu %>%
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

ebu_arhennius_FO <- function(A, a, temp, Q, P){
  est = A * a ^ (temp-20) + rnorm(100,0, sd = Q)
  return(est)
}

parms <- sample_n(res_coefficients, 100, replace=TRUE)
temp_data <- c(ebu_base_temp$temp_for_model_C)

# Process
out <- list()

for(s in 1:length(temp_data)){
  
  prediction <- ebu_arhennius_FO(temp = temp_data[s],
                                 A = mean(parms$A),
                                 a = mean(parms$a),
                                 Q = parms$sd.pro)
  out[[s]] <- prediction
}

prediction_output = as.data.frame(do.call(rbind, out))

prediction_output_long <- prediction_output %>% t(.) %>% reshape2::melt(.) %>%
  group_by(Var2) %>%
  summarize(mean = mean(value),
            sd = sd(value),
            var = var(value)) %>%
  mutate(mean = ifelse(mean <= 0, 0, mean))

var_pro = mean(prediction_output_long$sd)

# Process

out <- list()

for(s in 1:length(temp_data)){
  
  prediction <- ebu_arhennius_FO(temp = temp_data[s],
                                 A = parms$A,
                                 a = parms$a,
                                 Q = parms$sd.pro*0)
  out[[s]] <- prediction
}

prediction_output = as.data.frame(do.call(rbind, out))

prediction_output_long <- prediction_output %>% t(.) %>% reshape2::melt(.) %>%
  group_by(Var2) %>%
  summarize(mean = mean(value),
            sd = sd(value),
            var = var(value)) %>%
  mutate(mean = ifelse(mean <= 0, 0, mean))

var_par = mean(prediction_output_long$sd)

sum_var = var_par+var_pro

a_prop_par = var_par/sum_var
b_prop_pro = var_pro/sum_var

prop_var <- as.data.frame(rbind(a_prop_par,b_prop_pro,sum_var)) %>%
  mutate(type = "reservoir calibration") %>%
  mutate(waterbody_type = "reservoir") %>%
  rownames_to_column(., "row_names")


coefficients_all <- rbind(all_coefficients_DF, lake_coefficients_DF, res_coefficients_DF)

