s = Sys.time()

library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(vroom, warn.conflicts = FALSE)
library(readr, warn.conflicts = FALSE)


FOA_model <- function(temp){
  96.98125 * 1.17077 ^ (temp-20)
}

Aben_model <- function(temp){
  100 * 1.1 ^ (temp-20)
}

ebu_arhennius_FO <- function(A, a, temp, Q, P){
  est = A * a ^ (temp-20) + rnorm(100,0, sd = Q) + rnorm(100,0, sd = P)
  return(est)
}

file1 <- "/Users/ryanmcclure/Documents/GLEE-v1.0/data/GLEE-upscale-countries/Brazil.csv"

file2 <- "/Users/ryanmcclure/Documents/GLEE-v1.0/data/GLEE-upscale-countries/Germany.csv"

output <- "/central/groups/carnegie_poc/rmcclure/project-GLEE/global_ebullition_emission.csv"

d <- vroom::vroom(paste0(file2), col_names = F) %>%
    filter(X1 >= 2002) %>%
    mutate(lake_area_m2 = X7*1000000) %>%
    select(-X7) %>%
    na.omit(.) %>%
    mutate(ebullition_FO_mg_m2_d = ifelse(X9 > 50, 0, mapply(FOA_model, X8-273.15)),
           ebullition_Aben_mg_m2_d = ifelse(X9 > 50, 0, mapply(Aben_model, X8-273.15)),
           #ebullition_JAGS_mg_m2_d = ifelse(X9 > 75, 0, mapply(Bayes_arhennius_FO, 96.6, 1.17, X8-273.15, 232)),
           #ebullition_JAGS_mg_m2_d = ifelse(ebullition_JAGS_mg_m2_d<0,0,ebullition_JAGS_mg_m2_d),
           ebullition_LAC_FO = ebullition_FO_mg_m2_d*lake_area_m2*0.2,
           ebullition_LAC_Aben = ebullition_Aben_mg_m2_d*lake_area_m2*0.2,
           #ebullition_LAC_JAGS = ebullition_JAGS_mg_m2_d*lake_area_m2*0.000001*0.2*30,
           'date' = lubridate::make_date(year = X1, month = X2))

d2 <- d %>% select(ebullition_FO_mg_m2_d, ebullition_Aben_mg_m2_d,ebullition_LAC_FO,ebullition_LAC_Aben, date) %>%
  group_by(date) %>%
  summarize(mean_FO = mean(ebullition_FO_mg_m2_d),
            sd_FO = sd(ebullition_FO_mg_m2_d),
            mean_Aben = mean(ebullition_Aben_mg_m2_d),
            sd_Aben = sd(ebullition_Aben_mg_m2_d),
            mean_FO_LAC = mean(ebullition_LAC_FO),
            sd_FO_LAC = sd(ebullition_LAC_FO),
            mean_Aben_LAC = mean(ebullition_LAC_Aben),
            sd_Aben_LAC = sd(ebullition_LAC_Aben)) %>%
  mutate(sd_low_FO = mean_FO,
         sd_low_Aben = mean_Aben,
         sd_low_FO_LAC= mean_FO_LAC,
         sd_low_Aben_LAC = mean_Aben_LAC) %>%
  reshape2::melt(., id = c("date","sd_FO", "sd_Aben", "sd_low_FO", "sd_low_Aben", "sd_low_FO_LAC", "sd_low_Aben_LAC", "sd_FO_LAC", "sd_Aben_LAC"))

p3 <- d2 %>% filter(variable == "mean_FO_LAC") %>%
  mutate(value = value * 0.000001,
         sd_FO_LAC = sd_FO_LAC * 0.000001,
         sd_low_FO_LAC = sd_low_FO_LAC* 0.000001) %>%
  ggplot(., aes(x = date, y = value))+
  geom_line(lwd = 0.5, alpha = 1)+
  geom_ribbon(aes(x = date, y = value, ymin = value-sd_low_FO_LAC, ymax = value+sd_FO_LAC), alpha = 0.2, lwd = 0)+
  labs(title = "Assuming **NO** Observation or Process Error - Mean Values")+
  ylab("kg CH4/day")+
  ylim(0,350)+
  theme_classic()

p4 <- d2 %>% filter(variable == "mean_Aben_LAC") %>%
  mutate(value = value * 0.000001,
         sd_Aben_LAC = sd_Aben_LAC * 0.000001,
         sd_low_Aben_LAC = sd_low_Aben_LAC* 0.000001) %>%
  ggplot(., aes(x = date, y = value))+
  geom_line(lwd = 0.5, alpha = 1)+
  geom_ribbon(aes(x = date, y = value, ymin = value-sd_low_Aben_LAC, ymax = value+sd_Aben_LAC), alpha = 0.2, lwd = 0)+
  labs(title = "Assuming **NO** Observation or Process Error - Johnson Paper")+
  ylab("kg CH4/day")+
  ylim(0,350)+
  theme_classic()

d2 %>% filter(variable == "mean_FO_LAC") %>%
  summarise(mean_value = sum(value * 0.000001*30)/13.74,
            sd_FO_LAC = sum(sd_FO_LAC * 0.000001*30)/13.74,
            sd_low_FO_LAC = sum((value-sd_low_FO_LAC)* 0.000001*30)/13.74)

d2 %>% filter(variable == "mean_Aben_LAC") %>%
  summarise(mean_value = sum(value * 0.000001*30)/13.74,
            sd_Aben_LAC = sum(sd_Aben_LAC * 0.000001*30)/13.74,
            sd_low_Aben_LAC = sum((value-sd_low_Aben_LAC)* 0.000001*30)/13.74)





ggplot(d2, aes(x = date, y = value, group = variable, color = variable))+
  geom_line(lwd = 0.5, alpha = 1)+
  geom_point(data = ebu_obs_compare, aes(date, ch4_ebu), pch = 21, color = "grey", 
             fill = "red", inherit.aes = F)+
  labs(title = paste0("Number waterbodies plotted = ",length(unique(d$X3))))+
  theme_classic()


ebu_obs_compare <- ebu_base_temp %>%
  filter(country == "Germany") %>%
  mutate('date' = lubridate::make_date(year = year, month = month)) %>%
  filter(year > 2000)

length(unique(d$X3))

ggplot(d, aes(x = date, y = ebullition_JAGS_mg_m2_d, group = X3))+
  geom_line(lwd = 0.01, color = "green4", alpha = 0.5)+
   geom_point(data = ebu_obs_compare, aes(date, ch4_ebu), pch = 21, color = "grey", 
              fill = "red", inherit.aes = F)+
  labs(title = paste0("Number waterbodies plotted = ",length(unique(d$X3))))+
  theme_classic()



ebu_arhennius_FO <- function(A, a, temp, Q, P){
  est = A * a ^ (temp-20) + rnorm(100,0, sd = Q) + rnorm(100,0, sd = P)
  return(est)
}

parms <- sample_n(coefficients, 50, replace=TRUE)
temp_data <- c(ebu_base_temp$temp_for_model_C)

out <- list()

for(s in 1:length(temp_data)){
  
  prediction <- ebu_arhennius_FO(temp = temp_data[s],
                                 A = mean(parms$A),
                                 a = mean(parms$a),
                                 Q = parms$sd.pro*0,
                                 P = parms$sd.obs)
  out[[s]] <- prediction
}

prediction_output = as.data.frame(do.call(rbind, out))



e <- Sys.time()
t=e-s
print(t)


