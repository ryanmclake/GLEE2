


w <- f %>%
  select(ch4_diff, ch4_ebu, mean_annual_temp_k, mean_obs_wtemp_k, elevation, lat, waterbody_type, country, biome_type, hybas_id) %>%
  mutate(temp_for_model = ifelse(is.na(mean_annual_temp_k),mean_obs_wtemp_k, mean_annual_temp_k)) 


ebu <- w %>% select(country, biome_type, mean_ebu, ebu_sd, temp_for_model, waterbody_type, lat) %>% na.omit(.) %>% filter(mean_ebu >= 0) %>%
  mutate(temp_for_model = ifelse(temp_for_model < 200, NA, temp_for_model)) %>% 
  filter(waterbody_type != "river") %>%
  mutate(inverse_temp = 1/temp_for_model) %>%
  mutate(mean_ebu = log((mean_ebu))) %>%
  mutate(sd = 0.05) %>%
  na.omit(.) %>%
  filter(!is.infinite(mean_ebu))

ggplot(ebu, aes(x = inverse_temp, y = mean_ebu))+
  geom_point(aes(color = waterbody_type))+geom_smooth(method = "lm")+
  geom_abline(intercept = 17.59365, slope= -4087.339,linewidth=2)

library(ggplot2)

#define equation
my_equation <- function(x){(1.0/8.314)*(1/x) + 2.3}

#plot equation
ggplot(data.frame(x=c(1, 50)), aes(x=x)) + 
  stat_function(fun=my_equation)


model <- lm(mean_ebu~inverse_temp, data = ebu)
summary(model)




hist(ebu$mean_ebu, breaks = 100)

diff <- w %>% select(country, biome_type, mean_diff, diff_sd, temp_for_model, lat) %>% na.omit(.) %>% filter(mean_diff >= 0) %>%
  mutate(temp_for_model = ifelse(temp_for_model < -350, NA, temp_for_model)) %>% 
  mutate(mean_diff = log(mean_diff+0.1)) %>%
  mutate(sd = 0.05) %>%
  na.omit(.) %>%
  filter(mean_diff < 5000)

hist(diff$mean_diff, breaks = 100)


model= ("ebu_model.txt")
jagsscript = cat("
model {  
   
   #priors===================================================
   
   EA ~ dunif(0,1000)
   theta ~ dunif(0,100)
   sd.pro ~ dunif(0, 1000)
   
   #end priors===============================================
   
   for(i in 1:N) {
     
      #process model=============================================
     
      tau.pro[i] <- 1/((sd.pro)*(sd.pro))
      predX[i] <- (Ea/8.314)*(1/(temp[i]+A))
      X[i] ~ dnorm(predX[i],tau.pro[i])
     
      #end of process model======================================
     
      #data model================================================
     
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
            #end of data model=========================================
   }
  }", file = model)

# 
# model= ("ebu_model.txt")
# jagsscript = cat("
# model {  
#    
#    #priors===================================================
#    
#    ebtwenty ~ dunif(0,1000)
#    theta ~ dunif(0,100)
#    sd.pro ~ dunif(0, 1000)
#    
#    #end priors===============================================
#    
#    for(i in 1:N) {
#      
#       #process model=============================================
#      
#       tau.pro[i] <- 1/((sd.pro)*(sd.pro))
#       predX[i] <- ebtwenty * (theta^(temp[i]-20))
#       X[i] ~ dnorm(predX[i],tau.pro[i])
#      
#       #end of process model======================================
#      
#       #data model================================================
#      
#       Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
#             #end of data model=========================================
#    }
#   }", file = model)


set.seed(1)

# cal <- eb_mid1 %>% sample_n(., 100)
# cal_vals <- c(cal$row_names)
# val <- eb_mid1 %>% filter(!row_names %in% cal_vals) %>%
#   select(-row_names) %>%
#   tibble::rownames_to_column(., "row_names") %>%
#   mutate(row_names = as.numeric(row_names))


jags.data = list(Y = ebu$mean_ebu,
                 tau.obs = 1/((ebu$sd)) ^ 2,
                 N = nrow(ebu),
                 temp = ebu$temp_for_model)

nchain = 3
chain_seeds <- c(200,800,1400)
init <- list()
for(i in 1:nchain){
  init[[i]] <- list(sd.pro = runif(1, 0.01, 2),
                    Ea = runif(1, -3000,150),
                    A = runif(1, -100, 1000),
                    .RNG.name = "base::Wichmann-Hill",
                    .RNG.seed = chain_seeds[i])
}

j.model   <- jags.model(file = model,
                        data = jags.data,
                        inits = init,
                        n.chains = 3)

eval_ebu  <- coda.samples(model = j.model,
                          variable.names = c("sd.pro","Ea","A"),
                          n.iter = 200000, n.burnin = 20000, thin = 200)

plot(eval_ebu)
print(gelman.diag(eval_ebu))

parameter <- eval_ebu %>%
  spread_draws(sd.pro, Ea, A) %>%
  select(sd.pro, Ea, A)

parameter_new <- parameter %>%
  mutate(Ea = exp(Ea),
         A = exp(A))

exp(mean(parameter_new$Ea*8.314))

exp(mean(parameter_new$A))

hist(parameter$Ea, breaks = 100)
hist(parameter$A, breaks = 100)

summary(parameter)



temp_ebu_function_mcclure_mean <- function(x) {(-0.0008148/8.314*x) + exp(0.10124)}
temp_ebu_function_mcclure_min <- function(x) {(-0.000001/8.314*x) + exp(0.12357)}
temp_ebu_function_mcclure_max <- function(x) {(-0.0225249/8.314*x) + exp(0.14497)}

ebu_function_output <- data.frame(x = 253:313,            # Create data for ggplot2
                                  values = c(
                                             temp_ebu_function_mcclure_mean(253:313),
                                             temp_ebu_function_mcclure_min(253:313),
                                             temp_ebu_function_mcclure_max(253:313)),
                                  model = rep(c(
                                                "Bayes - Mean Parameter",
                                                "Bayes - 1st Q Parameters",
                                                "Bayes - 3rd Q Parameter"), each = 61))

ggplot(ebu, aes(x = 1/temp_for_model, y = mean_ebu))+
  geom_point()+
  geom_line(data = ebu_function_output, aes(1/x, values, group = model, color = model), lwd = 2)+
  scale_color_viridis(discrete = T, option = "C")+
  theme_classic()


ebu_validation <- function(ebtwenty, theta, temp, Q){
  est = ebtwenty * (theta ^ (temp - 20)) + rnorm(100,0, sd = Q)
  return(est)
}

parms <- sample_n(parameter, 1000, replace=TRUE)

out <- list()

for(s in 1:length(ebu$temp_for_model)){
  
  validation <- ebu_validation(temp = ebu$temp_for_model,
                               ebtwenty = parms$ebtwenty,
                               theta = parms$theta,
                               Q = mean(parms$sd.pro))
  out[[s]] <- validation
}

model_validate = as.data.frame(do.call(rbind, out)) %>%
  tibble::rownames_to_column(., "row_names") %>%
  pivot_longer(!row_names, names_to = "iteration", values_to = "value") %>%
  group_by(row_names) %>%
  summarize(mean_eb_flux = mean(value),
            sd_eb_flux = sd(value), 
            var_flux = var(value)) %>%
  arrange(row_names) %>%
  left_join(., ebu, by = "row_names")

mean <- as.data.frame(rowMeans(model_validate))

bind_cols(data_val$eb_flux, mean$`rowMeans(model_validate)`) %>%
  summarize(rmse = rmse(`...1`, `...2`))




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
                 tau.obs = 1/((diff$sd)) ^ 2,
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
temp_diff_function_mcclure_min <- function(x) {7.7140*exp(0.11034*x)}
temp_diff_function_mcclure_max <- function(x) {48.3526*exp(0.23625*x)}

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


