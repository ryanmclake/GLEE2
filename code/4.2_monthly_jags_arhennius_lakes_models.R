set.seed(643)

if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, MCMCvis, lubridate, tidybayes,
               ncdf4, reshape2, zoo, patchwork, hydroGOF, viridis,
               imputeTS, devtools, scales, forecast, coda, rjags, R2jags,gridExtra)

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
      predX[i] <- A * a ^ (temp[i] - 20)
      X[i] ~ dnorm(predX[i],tau.pro[i])
     
      #end of process model======================================
     
      #data model================================================
     
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
            #end of data model=========================================
   }
  }", file = ahrennius_ebu_model)

ebu_coefficients <- read_csv("./data/ebullition_coefficients.csv")

summary(ebu_coefficients)

### MONTHLY TIMESTEP ###
ebu_base_temp <- base %>% select(ch4_ebu, temp_for_model_K, waterbody_type, waterbody_id, month, tot_sampling_events) %>%
  na.omit(.) %>%
  filter(waterbody_type == "lake") %>%
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
                    A = runif(1, 728.1, sd(ebu_coefficients$E20)),
                    a = runif(1, 1.136, 3),
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


coefficients <- eval_ebu %>%
  spread_draws(sd.pro, A, a) %>%
  filter(.chain == 1)

summary(coefficients)

ebu_arhennius_FO <- function(A, a, temp, Q){
  est = A * a ^ (temp-20) + rnorm(100,0, sd = Q)
  return(est)
}

parms <- sample_n(coefficients, 100, replace=TRUE)

mean(ebu_arhennius_FO(temp = 30,
                      A = parms$A,
                      a = parms$a,
                      Q = parms$sd.pro))

temp_data <- c(ebu_base_temp$temp_for_model_C)

out <- list()

for(s in 1:length(temp_data)){
  
  prediction <- ebu_arhennius_FO(temp = temp_data[s],
                                 A = parms$A,
                                 a = parms$a,
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

validate <- cbind(prediction_output_long, ebu_base_temp$temp_for_model_C, ebu_base_temp$ch4_ebu) %>%
  summarise(NSE = NSE(mean, `ebu_base_temp$ch4_ebu`),
            rmse = rmse(mean, `ebu_base_temp$ch4_ebu`))

#var_E20 = mean(prediction_output_long$var)

#var_theta = mean(prediction_output_long$var)

#var_process = mean(prediction_output_long$var)

sum_var = var_E20+var_theta+var_process

a_prop_E20 = var_E20/sum_var
b_prop_theta = var_theta/sum_var
c_prop_process = var_process/sum_var

prop_var <- as.data.frame(rbind(a_prop_E20,b_prop_theta,c_prop_process)) %>%
  mutate(type = "monthly calibration") %>%
  rownames_to_column(., "row_names")

ggplot(prop_var, aes(x = type, y = V1, fill = row_names)) +
  geom_bar(stat = 'identity', position = 'stack')


