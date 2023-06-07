
# Clean your environment
rm(list=ls())

#Clean your memory
gc()

# Set seed for repeatability
set.seed(3818)

# install packages needed for the analysis

if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, MCMCvis, lubridate, tidybayes, readr,
               ncdf4, reshape2, zoo, patchwork, hydroGOF, viridis,
               imputeTS, devtools, scales, forecast, coda, rjags, 
               R2jags, gridExtra, maps, hexbin, rnaturalearth, sf)

# Load in the prior Ebullition Estimates
ebu_coefficients <- read_csv("./source_data/ebullition_coefficients.csv")

# obtain the summary statistics around our ebullition priors
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

# Load in the data base with the assigned error
error_ebu <- readr::read_csv("./output/filtered_GLEE_ebullition_w_GLCP_link_and_ERROR.csv") %>%
  filter(waterbody_type == "lake" | waterbody_type == "reservoir") %>%
  filter(sd_time > 0)

# specify the lake type that we will filter to calibrate
# we will calibrate lakes and reservoirs independently
lake_type <- c(unique(error_ebu$waterbody_type))

# specify our data frames that we want out
out_raw_coefficients <- list()
out_agg_coefficients <- list()

# for loop that loops between lakes and reservoirs
for(g in 1:length(lake_type)){
  
  # get either a lake or reservoir
  error_ebu_new <- error_ebu %>% 
    filter(waterbody_type == lake_type[g])
  
  # specify our data for JAGS
  jags.data = list(Y = error_ebu_new$ch4_ebu,
                   #tau.obs = 1/(error_ebu_new$sd_time) ^ 2,
                   #tau.obs = 1/(error_ebu_new$sd_space_time) ^ 2,
                   tau.obs = 1/(error_ebu_new$sd_space) ^ 2,
                   N = nrow(error_ebu_new),
                   temp = as.numeric(error_ebu_new$temp_for_model_C))
  
  # Set up our chains
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
  
  
  # prepare the model
  j.model   <- jags.model(file = arrhenius_ebu_model,
                          data = jags.data,
                          inits = init,
                          n.chains = 3)
  # evaluate the MCMCs
  eval_ebu  <- coda.samples(model = j.model,
                            variable.names = c("sd.pro","E20","omega"),
                            n.iter = 20000, n.burnin = 2000, thin = 20)
  
  #Plot the output of the MCMC
  plot(eval_ebu)
  # obtain the gelman diagnostics
  print(gelman.diag(eval_ebu))
  
  #save the raw coefficient values form the MCMC run
  coefficients_raw <- eval_ebu %>%
    spread_draws(sd.pro, E20, omega) %>%
    mutate(waterbody_type = lake_type[g])
  
  # I do this because it is a list object before the for loop
  out_raw_coefficients[[g]] <- coefficients_raw
    
  # this aggregates the coefficeients to get the mean, SD and 95% CIs of the coefficients
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

# output the raw and aggregated coefficients as a DF
out_raw_coefficients = as.data.frame(do.call(rbind, out_raw_coefficients))
out_agg_coefficients = as.data.frame(do.call(rbind, out_agg_coefficients))

# generate a posterior prediction function
arrhenius_ebu_posterior_prediction <- function(E20, omega, temp, Q){
  est = (E20 * omega ^ (temp-20)) + rnorm(1000, 0, sd = Q)
  return(est)
}

# build to empty data frames to fill in for the model validation.
validate_models <- list()
validation <- list()

# another loop that goes between Lakes and Reservoirs 
for(h in 1:length(lake_type)){

  # randomly sample 1000 coefficient values
parms <- out_raw_coefficients %>%
  filter(waterbody_type == lake_type[h]) %>%
  sample_n(., 1000, replace=TRUE)

  # pull the temperature data from the original DF 
temp_data <- error_ebu %>% 
  filter(waterbody_type == lake_type[h]) %>%
  select(temp_for_model_C)

  # concatinate the temp data
temp_data <- c(temp_data$temp_for_model_C)

  # pull the original observed data to compare against the predictions
ebu_data <- error_ebu %>% 
  filter(waterbody_type == lake_type[h]) %>%
  select(ch4_ebu, temp_for_model_C, waterbody_type)

  # build a loop that goes through all 1000 coefficient iterations I randomly pulled
predict_output <- list()

for(s in 1:length(temp_data)){
  
  prediction <- arrhenius_ebu_posterior_prediction(temp = temp_data[s],
                                                   E20 = parms$E20,
                                                   omega = parms$omega,
                                                   Q = parms$sd.pro)
  predict_output[[s]] <- prediction
}

# output the predictions as a DF
prediction_output = as.data.frame(do.call(rbind, predict_output))

# aggregate the prediction output and summarize the predictions by mean, SD , and var
prediction_output_long <- prediction_output %>% t(.) %>% reshape2::melt(.) %>%
  group_by(Var2) %>%
  summarize(mean = mean(value),
            sd = sd(value),
            var = var(value))

# see how well the model does against the obervations
validate <- cbind(prediction_output_long, ebu_data)

validation[[h]] <- validate

# generate an old school plot of the model output and the observations
# the plot is basically plotting the equation 
# where temp is on the X and ebullition observed and predcited is on the Y
plot(validate$temp_for_model_C, validate$ch4_ebu, col = "blue4", pch = 16)
points(validate$temp_for_model_C, validate$mean, col = "black", pch = 21, bg = "red")

# determined the NSE and RMSE of the emissions in lakes and reservoirs
NSE_RMSE <- cbind(prediction_output_long, ebu_data) %>%
  summarise(NSE = NSE(mean, ch4_ebu),
            rmse = rmse(mean, ch4_ebu))

validate_models[[h]] <- NSE_RMSE

}

# Get the model evaluation (NSE & RMSE) output
validate_models = as.data.frame(do.call(rbind, validate_models))
validation = as.data.frame(do.call(rbind, validation))


