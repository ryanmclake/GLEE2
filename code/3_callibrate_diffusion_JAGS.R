# Clean your environment
rm(list=ls())

#Clean your memory
gc()

# Set seed for repeat measurements
set.seed(3818)

# install packages needed for the analysis

if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, MCMCvis, lubridate, tidybayes, readr,
               ncdf4, reshape2, zoo, patchwork, hydroGOF, viridis,
               imputeTS, devtools, scales, forecast, coda, rjags, 
               R2jags, gridExtra, maps, hexbin, rnaturalearth, sf)


#### NOTE ######
# EVERYTHING BELOW IS THE SAME AS THE 2_calibrate_ebullition_JAGS.R script except the 
# model structure in the JAGS model. Consequently, I have not annotated it as well because 
# I honestly did not want to retype everything twice. 

# Load in the prior diffusion Estimates
diff_coefficients <- read_csv("./source_data/diffusion_coefficients.csv")

summary(diff_coefficients)
sd(diff_coefficients$a)
sd(diff_coefficients$B)

# Build the JAGS MODEL
arrhenius_diff_model = ("./output/arrhenius_diff_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   A ~ dnorm(0.07350,1/0.07288278)
   a ~ dnorm(0.12683,1/0.043545)
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
  }", file = arrhenius_diff_model)


error_diff <- readr::read_csv("./output/filtered_GLEE_diffusion_w_GLCP_link_and_ERROR.csv") %>%
  filter(waterbody_type == "lake" | waterbody_type == "reservoir") %>%
  filter(sd_time > 0)

lake_type <- c(unique(error_diff$waterbody_type))

out_raw_coefficients <- list()
out_agg_coefficients <- list()

for(g in 1:length(lake_type)){
  
  error_diff_new <- error_diff %>% 
    filter(waterbody_type == lake_type[g])
  
  jags.data = list(Y = error_diff_new$ch4_diff,
                   #tau.obs = 1/(error_diff_new$sd_time) ^ 2,
                   #tau.obs = 1/(error_diff_new$sd_space_time) ^ 2,
                   tau.obs = 1/(error_diff_new$sd_space) ^ 2,
                   N = nrow(error_diff_new),
                   temp = as.numeric(error_diff_new$temp_for_model_C))
  
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(sd.pro = runif(1, 1, 3),
                      A = runif(1, 0.07350, 2),
                      a = runif(1, 0.12683, 2),
                      .RNG.name = "base::Wichmann-Hill",
                      .RNG.seed = chain_seeds[i])
  }
  
  
  
  j.model   <- jags.model(file = arrhenius_diff_model,
                          data = jags.data,
                          inits = init,
                          n.chains = 3)
  
  eval_diff  <- coda.samples(model = j.model,
                            variable.names = c("sd.pro","A","a"),
                            n.iter = 20000, n.burnin = 2000, thin = 20)
  
  plot(eval_diff)
  print(gelman.diag(eval_diff))
  
  coefficients_raw <- eval_diff %>%
    spread_draws(sd.pro, A, a) %>%
    mutate(waterbody_type = lake_type[g])
  
  out_raw_coefficients[[g]] <- coefficients_raw
  
  
  coefficients_aggregate <- coefficients_raw %>%
    summarise(mean_A = mean(A),
              sd_A = sd(A),
              upper_95_A = quantile(A, 0.95, na.rm = T),
              lower_95_A = quantile(A, 0.05, na.rm = T),
              mean_a = mean(a),
              sd_a = sd(a),
              upper_95_a = quantile(a, 0.95, na.rm = T),
              lower_95_a = quantile(a, 0.05, na.rm = T),
              mean_process = mean(sd.pro),
              sd_process = sd(sd.pro),
              upper_95_process = quantile(sd.pro, 0.95, na.rm = T),
              lower_95_process = quantile(sd.pro, 0.05, na.rm = T)) %>%
    mutate(waterbody_type = lake_type[g])
  
  out_agg_coefficients[[g]] <- coefficients_aggregate
  
}


out_raw_coefficients = as.data.frame(do.call(rbind, out_raw_coefficients))
write_csv(out_raw_coefficients, "./output/raw_diffusion_coefficients.csv")

out_agg_coefficients = as.data.frame(do.call(rbind, out_agg_coefficients))


arrhenius_diff_posterior_prediction <- function(A, a, temp, Q){
  est = (A*exp(a*(temp))) + rnorm(1000, 0, sd = Q)
  return(est)
}

validate_models <- list()
validation <- list()

for(h in 1:length(lake_type)){
  
  parms <- out_raw_coefficients %>%
    filter(waterbody_type == lake_type[h]) %>%
    sample_n(., 1000, replace=TRUE)
  
  temp_data <- error_diff %>% 
    filter(waterbody_type == lake_type[h]) %>%
    select(temp_for_model_C)
  
  temp_data <- c(temp_data$temp_for_model_C)
  
  diff_data <- error_diff %>% 
    filter(waterbody_type == lake_type[h]) %>%
    select(ch4_diff, temp_for_model_C, waterbody_type)
  
  predict_output <- list()
  
  for(s in 1:length(temp_data)){
    
    prediction <- arrhenius_diff_posterior_prediction(temp = temp_data[s],
                                                     A = parms$A,
                                                     a = parms$a,
                                                     Q = parms$sd.pro)
    predict_output[[s]] <- prediction
  }
  
  prediction_output = as.data.frame(do.call(rbind, predict_output))
  
  prediction_output_long <- prediction_output %>% t(.) %>% reshape2::melt(.) %>%
    group_by(Var2) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              var = var(value))
  
  validate <- cbind(prediction_output_long, diff_data)
  
  validation[[h]] <- validate
  
  plot(validate$temp_for_model_C, validate$ch4_diff, col = "blue4", pch = 16)
  points(validate$temp_for_model_C, validate$mean, col = "black", pch = 21, bg = "red")
  
  NSE_RMSE <- cbind(prediction_output_long, diff_data) %>%
    summarise(NSE = NSE(mean, ch4_diff),
              rmse = rmse(mean, ch4_diff))
  
  validate_models[[h]] <- NSE_RMSE
  
}

validate_models = as.data.frame(do.call(rbind, validate_models))
validation = as.data.frame(do.call(rbind, validation))
