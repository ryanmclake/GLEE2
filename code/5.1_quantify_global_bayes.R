library(tidyverse)

parameters

ebu_arhennius_FO_1_4 <- function(temp){
  est = rnorm(1,95.89238, sd = 24.677516) * rnorm(1,1.175983, sd = 0.05976376) ^ (temp-20) + rnorm(1,0, sd = 56.98130) + rnorm(1,0, sd = 224.1444660)
  return(est)
}

ebu_arhennius_FO_50 <- function(temp){
  est = rnorm(1,97.12485, sd = 8.341771) * rnorm(1,1.171232, sd = 0.01409348) ^ (temp-20) + rnorm(1,0, sd = 99.93237) + rnorm(1,0, sd = 0.9928644)
  return(est)
}

out_jags <- list()
iteration <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20)

for(i in 1:length(iteration)){
  
  d <- vroom::vroom(paste0(file2), col_names = F) %>%
    filter(X1 >= 2002) %>%
    mutate(lake_area_m2 = X7*1000000) %>%
    select(-X7) %>%
    na.omit(.) %>%
    mutate(ebu_JAGS_FO_1_4 = ifelse(X9 > 50, 0, mapply(ebu_arhennius_FO_1_4, X8-273.15)),
           ebu_JAGS_FO_50 = ifelse(X9 > 50, 0, mapply(ebu_arhennius_FO_50, X8-273.15)),
           ebullition_LAC_1_4 = ebu_JAGS_FO_1_4*lake_area_m2*0.4,
           ebullition_LAC_50 = ebu_JAGS_FO_50*lake_area_m2*0.4,
           'date' = lubridate::make_date(year = X1, month = X2),
           iteration = iteration[i])
  
  out_jags[[i]] <- d
  
}

JAGS_output = as.data.frame(do.call(rbind, out_jags))

JAGS_output_summary <- JAGS_output %>%
  group_by(date) %>%
  summarize(mean_1_4 = mean(ebu_JAGS_FO_1_4),
            sd_1_4 = sd(ebu_JAGS_FO_1_4),
            mean_50 = mean(ebu_JAGS_FO_50),
            sd_50 = sd(ebu_JAGS_FO_50),
            mean_1_4_LAC = mean(ebullition_LAC_1_4),
            sd_1_4_LAC = sd(ebullition_LAC_1_4),
            mean_50_LAC = mean(ebullition_LAC_50),
            sd_50_LAC = sd(ebullition_LAC_50)) %>%
  mutate(sd_low_1_4 = mean_1_4,
         sd_low_50 = mean_50,
         sd_low_1_4_LAC = mean_1_4_LAC,
         sd_low_50_LAC = mean_50_LAC) %>%
  reshape2::melt(., id = c("date", "sd_50", "sd_1_4", "sd_low_1_4", "sd_low_50", "sd_50_LAC", "sd_1_4_LAC", "sd_low_1_4_LAC", "sd_low_50_LAC"))

JAGS_output_summary %>% filter(variable == "mean_1_4") %>%
ggplot(., aes(x = date, y = value))+
  geom_line(lwd = 0.5, alpha = 1)+
  geom_ribbon(aes(x = date, y = value, ymin = value-sd_low_1_4, ymax = value+sd_1_4), alpha = 0.2, lwd = 0)+
  geom_point(data = ebu_obs_compare, aes(date, ch4_ebu), pch = 21, color = "grey", 
             fill = "red", inherit.aes = F)+
  labs(title = "Assuming Low Observation Certainty (i.e. 1-4 total samples observed)")+
  theme_classic()

p1 <- JAGS_output_summary %>% filter(variable == "mean_1_4_LAC") %>%
  mutate(value = ifelse(value < 0, 0, value * 0.000001),
         sd_1_4_LAC = sd_1_4_LAC * 0.000001,
         sd_low_1_4_LAC = sd_low_1_4_LAC* 0.000001) %>%
  ggplot(., aes(x = date, y = value))+
  geom_line(lwd = 0.5, alpha = 1)+
  geom_ribbon(aes(x = date, y = value, ymin = value-sd_low_1_4_LAC, ymax = value+sd_1_4_LAC), alpha = 0.2, lwd = 0)+
  labs(title = "Includes Low Observation Certainty (1-4) and Model Process Error")+
  ylab("kg CH4/day")+
  ylim(0,350)+
  theme_classic()

JAGS_output_summary %>% filter(variable == "mean_50") %>%
  ggplot(., aes(x = date, y = value))+
  geom_line(lwd = 0.5, alpha = 1)+
  geom_ribbon(aes(x = date, y = value, ymin = value-sd_low_50, ymax = value+sd_50), alpha = 0.2, lwd = 0)+
  geom_point(data = ebu_obs_compare, aes(date, ch4_ebu), pch = 21, color = "grey", 
             fill = "red", inherit.aes = F)+
  labs(title = "Assuming High Observation Certainty (i.e. 50+ samples) and Includes Model Process Error")+
  ylim(0,350)+
  theme_classic()

p2 <- JAGS_output_summary %>% filter(variable == "mean_50_LAC") %>%
  mutate(value = ifelse(value < 0, 0, value * 0.000001),
         sd_50_LAC = sd_50_LAC * 0.000001,
         sd_low_50_LAC = sd_low_50_LAC* 0.000001) %>%
  ggplot(., aes(x = date, y = value))+
  geom_line(lwd = 0.5, alpha = 1)+
  geom_ribbon(aes(x = date, y = value, ymin = value-sd_low_50_LAC, ymax = value+sd_50_LAC), alpha = 0.2, lwd = 0)+
  labs(title = "Includes High Observation Certainty (50+) and Model Process Error")+
  ylab("kg CH4/day")+
  ylim(0,350)+
  theme_classic()

library(patchwork)

plots <- (p4+p1)/(plot_spacer()+p2)







JAGS_output_summary %>% filter(variable == "mean_50_LAC") %>%
  summarise(mean_value = sum(value * 0.000001*30)/13.74,
         sd_50_LAC_high = sum(sd_50_LAC * 0.000001*30)/13.74,
         sd_low_50_LAC = sum((value-sd_low_50_LAC)* 0.000001*30)/13.74)

JAGS_output_summary %>% filter(variable == "mean_1_4_LAC") %>%
  summarise(mean_value = sum(value * 0.000001*30)/13.74,
            sd_1_4_LAC_high = sum(sd_1_4_LAC * 0.000001*30)/13.74,
            sd_low_1_4_LAC = sum((value-sd_low_1_4_LAC)* 0.000001*30)/13.74)



