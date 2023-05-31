set.seed(643)

if (!"pacman" %in% installed.packages()) install.packages("pacman")
pacman::p_load(tidyverse, MCMCvis, lubridate, tidybayes,
               ncdf4, reshape2, zoo, patchwork, hydroGOF, viridis,
               imputeTS, devtools, scales, forecast, coda, rjags, R2jags,gridExtra)


sd_values = c(unique(diff_error2$se))

out_big <- list()
out_coefficients <- list()

for(k in 1:length(sd_values)){
  
  ahrennius_diff_model = ("arhennius_diff_model.txt")
  jagsscript = cat("
model {  
   
   #priors===================================================
   
   A ~ dnorm(0.07350,1/0.07288278)
   a ~ dnorm(0.12683,1/0.043545)
   sd.pro ~ dunif(0, 1000)
   sd.obs ~ dgamma(1/sd_values,1/sd_values)
   
   #end priors===============================================
   
   for(i in 1:N) {
     
      #process model=============================================
     
      tau.pro[i] <- 1/((sd.pro)*(sd.pro))
      predX[i] <- A*exp(a*(temp[i]))
      X[i] ~ dnorm(predX[i],tau.pro[i])
     
      #end of process model======================================
     
      #data model================================================
     
      Y[i] ~ dnorm(X[i], tau.obs[i]) # Observation variation
      tau.obs[i] <- 1/((sd.obs)*(sd.obs))
            #end of data model=========================================
   }
  }", file = ahrennius_diff_model)
  
  ### MONTHLY TIMESTEP ###
  diff_base_temp <- base %>% select(ch4_diff, temp_for_model_K,country, waterbody_type, year, month, tot_sampling_events) %>%
    na.omit(.) %>%
    mutate(temp_for_model_C = temp_for_model_K-273.15)
  
  jags.data = list(Y = diff_base_temp$ch4_diff,
                   N = nrow(diff_base_temp),
                   temp = as.numeric(diff_base_temp$temp_for_model_C),
                   sd_values = sd_values[k])
  
  nchain = 3
  chain_seeds <- c(200,800,1400)
  init <- list()
  for(i in 1:nchain){
    init[[i]] <- list(sd.pro = runif(1, 1, 3),
                      sd.obs = dgamma(1/sd_values[k], 1/sd_values[k]),
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
                            variable.names = c("sd.pro","A","a", "sd.obs"),
                            n.iter = 20000, n.burnin = 2000, thin = 20)
  plot(eval_diff)
  print(gelman.diag(eval_diff))
  
  coefficients <- eval_diff %>%
    spread_draws(sd.pro, A, a, sd.obs) %>%
    filter(.chain == 1)
  
  coefficients_DF <- coefficients %>%
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
              lower_95_process = quantile(sd.pro, 0.05, na.rm = T),
              mean_obs = mean(sd.obs),
              sd_obs = sd(sd.obs),
              upper_95_obs = quantile(sd.obs, 0.95, na.rm = T),
              lower_95_obs = quantile(sd.obs, 0.05, na.rm = T))%>%
    mutate(prediction_level = sd_values[k])
  
  out_coefficients[[k]] <- coefficients_DF
  
  
  diff_arhennius_FO <- function(A, a, temp, Q, P){
    est = A * exp(a*temp) + rnorm(100,0, sd = Q) + rnorm(100,0, sd = P)
    return(est)
  }
  
  parms <- sample_n(coefficients, 100, replace=TRUE)
  temp_data <- c(diff_base_temp$temp_for_model_C)
  
  out <- list()
  
  for(s in 1:length(temp_data)){
    
    prediction <- diff_arhennius_FO(temp = temp_data[s],
                                   A = mean(parms$A),
                                   a = mean(parms$a),
                                   Q = parms$sd.pro*0,
                                   P = parms$sd.obs)
    out[[s]] <- prediction
  }
  
  prediction_output = as.data.frame(do.call(rbind, out))
  
  prediction_output_long <- prediction_output %>% t(.) %>% reshape2::melt(.) %>%
    group_by(Var2) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              var = var(value)) %>%
    mutate(mean = ifelse(mean <= 0, 0, mean))
  
  var_obs = mean(prediction_output_long$sd)
  
  # Process
  
  out <- list()
  
  for(s in 1:length(temp_data)){
    
    prediction <- diff_arhennius_FO(temp = temp_data[s],
                                   A = mean(parms$A),
                                   a = mean(parms$a),
                                   Q = parms$sd.pro,
                                   P = parms$sd.obs*0)
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
  
  
  # Parameter
  
  out <- list()
  
  for(s in 1:length(temp_data)){
    
    prediction <- diff_arhennius_FO(temp = temp_data[s],
                                   A = parms$A,
                                   a = parms$a,
                                   Q = parms$sd.pro*0,
                                   P = parms$sd.obs*0)
    out[[s]] <- prediction
  }
  
  prediction_output = as.data.frame(do.call(rbind, out))
  
  prediction_output_long <- prediction_output %>% t(.) %>% reshape2::melt(.) %>%
    group_by(Var2) %>%
    summarize(mean = mean(value),
              sd = sd(value),
              var = var(value)) %>%
    mutate(mean = ifelse(mean <= 0, 0, mean))
  
  var_para = mean(prediction_output_long$sd)
  
  sum_var = var_para+var_pro+var_obs
  
  a_prop_para = var_para/sum_var
  b_prop_obs = var_obs/sum_var
  c_prop_pro = var_pro/sum_var
  
  prop_var <- as.data.frame(rbind(a_prop_para,b_prop_obs,c_prop_pro,sum_var)) %>%
    mutate(type = "monthly calibration") %>%
    mutate(se = sd_values[k]) %>%
    rownames_to_column(., "row_names")
  
  out_big[[k]] <- prop_var
  
  p <- ggplot(prop_var, aes(x = se, y = V1, fill = row_names)) +
    xlab(sd_values[k])+
    geom_bar(stat = 'identity', position = 'stack')
  
  ggsave(p, path = ".", filename = paste0("./uncertatinty_partition_",sd_values[k],"_.jpg"),
         width = 8, height = 4, device='jpg', dpi=175)
  
}

proportion_output = as.data.frame(do.call(rbind, out_big))
parameters_output = as.data.frame(do.call(rbind, out_coefficients))

proportions <- proportion_output %>% rename(se = se) %>% left_join(., diff_error2, by = "se") %>%
  filter(row_names != "sum_var")

parameters <- parameters_output %>% rename(se = prediction_level) %>% left_join(., diff_error2, by = "se")

coefficient <- proportion_output %>% rename(se = se) %>% left_join(., diff_error2, by = "se") %>%
  filter(row_names != "sum_var")

total_var <- proportion_output %>% rename(se = se) %>% left_join(., diff_error2, by = "se") %>%
  filter(row_names == "sum_var")

a <- proportions%>%
  group_by(row_names)%>% 
  ggplot(., aes(x = sample_size, y = V1, group=row_names, fill = row_names))+
  geom_area(aes(fill = row_names))+
  scale_fill_manual(values = c("#D55E00", "#CC79A7", "#56B4E9"))+
  theme_bw()+
  labs(title = "A: Partitioned uncertainty with increasing sample size")+
  ylab("Proportion of total variance")+
  xlab("")+
  theme(axis.text=element_text(size=15, color = "black"),
        axis.title=element_text(size=15, color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        title = element_text(size = 15), legend.position = "top",
        legend.text = element_text(size = 10, color = "black"))

b <- ggplot(total_var) +
  geom_point(aes(sample_size, V1), size = 10, pch = 21, color = "black", fill = "blue3")+
  ylab("Total Posterior Prediction SD")+
  xlab("Total samples from sites")+
  theme_classic()


