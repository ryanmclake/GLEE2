rm(list=ls())
gc()
set.seed(3818)

# install packages needed for the analysis

if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, MCMCvis, lubridate, tidybayes, readr,
               ncdf4, reshape2, zoo, patchwork, hydroGOF, viridis,
               imputeTS, devtools, scales, forecast, coda, rjags, 
               R2jags, gridExtra, maps, hexbin, rnaturalearth, sf)

# Load in the prior Ebullition Estimates
ebu_coefficients <- read_csv("./source_data/ebullition_coefficients.csv")

summary(ebu_coefficients)
sd(ebu_coefficients$omega)
sd(ebu_coefficients$E20)

# Build the JAGS MODEL
arrhenius_ebu_model = ("./output/arrhenius_ebu_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   E20 ~ dnorm(728.1,1/829.4559^2)
   omega ~ dnorm(1.136,1/0.128368^2)
   sd.pro ~ dunif(0,1000)
   
   #end priors===============================================
   
   for(i in 1:N) {
     
      #process model=============================================
     
      tau.pro[i] <- 1/((sd.pro)*(sd.pro))
      predX[i] <- E20 * (omega^(temp[i]-20))
      X[i] ~ dnorm(predX[i],tau.pro[i])
     
      #end of process model======================================
     
      #data model================================================
     
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
            #end of data model=========================================
   }
  }", file = arrhenius_ebu_model)


error_ebu <- readr::read_csv("./output/filtered_GLEE_ebullition_w_GLCP_link_and_ERROR.csv") %>%
  filter(waterbody_type == "lake" | waterbody_type == "reservoir") %>%
  filter(sd_time > 0)

lake_type <- c(unique(error_ebu$waterbody_type))

out_raw_coefficients <- list()
out_agg_coefficients <- list()

for(g in 1:length(lake_type)){
  
  error_ebu_new <- error_ebu %>% 
    filter(waterbody_type == lake_type[g])
  
  jags.data = list(Y = error_ebu_new$ch4_ebu,
                   tau.obs = 1/(error_ebu_new$sd_time) ^ 2,
                   N = nrow(error_ebu_new),
                   temp = as.numeric(error_ebu_new$temp_for_model_C))
  
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(sd.pro = runif(1, 1, 3),
                      E20 = runif(1, mean(ebu_coefficients$E20), sd(ebu_coefficients$E20)),
                      omega = runif(1, mean(ebu_coefficients$omega), 3),
                      .RNG.name = "base::Wichmann-Hill",
                      .RNG.seed = chain_seeds[i])
  }
  
  
  
  j.model   <- jags.model(file = arrhenius_ebu_model,
                          data = jags.data,
                          inits = init,
                          n.chains = 3)
  
  eval_ebu  <- coda.samples(model = j.model,
                            variable.names = c("sd.pro","E20","omega"),
                            n.iter = 20000, n.burnin = 2000, thin = 20)
  
  plot(eval_ebu)
  print(gelman.diag(eval_ebu))
  
  coefficients_raw <- eval_ebu %>%
    spread_draws(sd.pro, E20, omega) %>%
    mutate(waterbody_type = lake_type[g])
  
  out_raw_coefficients[[g]] <- coefficients_raw
    
  
  coefficients_aggregate <- coefficients_raw %>%
    summarise(mean_E20 = mean(E20),
              sd_E20 = sd(E20),
              upper_95_E20 = quantile(E20, 0.95, na.rm = T),
              lower_95_E20 = quantile(E20, 0.05, na.rm = T),
              mean_omega = mean(omega),
              sd_omega = sd(omega),
              upper_95_omega = quantile(omega, 0.95, na.rm = T),
              lower_95_omega = quantile(omega, 0.05, na.rm = T),
              mean_process = mean(sd.pro),
              sd_process = sd(sd.pro),
              upper_95_process = quantile(sd.pro, 0.95, na.rm = T),
              lower_95_process = quantile(sd.pro, 0.05, na.rm = T)) %>%
    mutate(waterbody_type = lake_type[g])
  
  out_agg_coefficients[[g]] <- coefficients_aggregate
  
}


out_raw_coefficients = as.data.frame(do.call(rbind, out_raw_coefficients))
out_agg_coefficients = as.data.frame(do.call(rbind, out_agg_coefficients))


arrhenius_ebu_posterior_prediction <- function(E20, omega, temp, Q){
  est = (E20 * omega ^ (temp-20)) + rnorm(1000, 0, sd = Q)
  return(est)
}

validate_models <- list()
validation <- list()

for(h in 1:length(lake_type)){

parms <- out_raw_coefficients %>%
  filter(waterbody_type == lake_type[h]) %>%
  sample_n(., 1000, replace=TRUE)

temp_data <- error_ebu %>% 
  filter(waterbody_type == lake_type[h]) %>%
  select(temp_for_model_C)

temp_data <- c(temp_data$temp_for_model_C)

ebu_data <- error_ebu %>% 
  filter(waterbody_type == lake_type[h]) %>%
  select(ch4_ebu, temp_for_model_C, waterbody_type)

predict_output <- list()

for(s in 1:length(temp_data)){
  
  prediction <- arrhenius_ebu_posterior_prediction(temp = temp_data[s],
                                                   E20 = parms$E20,
                                                   omega = parms$omega,
                                                   Q = parms$sd.pro)
  predict_output[[s]] <- prediction
}

prediction_output = as.data.frame(do.call(rbind, predict_output))

prediction_output_long <- prediction_output %>% t(.) %>% reshape2::melt(.) %>%
  group_by(Var2) %>%
  summarize(mean = mean(value),
            sd = sd(value),
            var = var(value))

validate <- cbind(prediction_output_long, ebu_data)

validation[[h]] <- validate

plot(validate$temp_for_model_C, validate$ch4_ebu, col = "blue4", pch = 16)
points(validate$temp_for_model_C, validate$mean, col = "black", pch = 21, bg = "red")

NSE_RMSE <- cbind(prediction_output_long, ebu_data) %>%
  summarise(NSE = NSE(mean, ch4_ebu),
            rmse = rmse(mean, ch4_ebu))

validate_models[[h]] <- NSE_RMSE

}

validate_models = as.data.frame(do.call(rbind, validate_models))
validation = as.data.frame(do.call(rbind, validation))
