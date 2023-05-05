
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(vroom, warn.conflicts = FALSE)
library(foreach, warn.conflicts = FALSE)
library(iterators, warn.conflicts = FALSE)
library(parallel, warn.conflicts = FALSE)
library(doParallel, warn.conflicts = FALSE)

s = Sys.time()

FOA_model <- function(temp){
  96.98125 * 1.17077 ^ (temp-20)
}

Aben_model <- function(temp){
  100 * 1.1 ^ (temp-20)
}

country <- list.files(path = "./data/GLEE-upscale-countries")
country <- gsub("\\..*", "", country)

emission_function <- function(x){
  
vroom::vroom(paste0("./data/GLEE-upscale-countries/",x,".csv"), col_names = F) %>%
    rename(year = X1, month = X2, hylak_id = X3, lat = X4, lon = X5, 
           country = X6, lake_area_km2 = X7, mean_temp_k = X8, percent_ice_cover = X9) %>%
    filter(year >= 2000) %>%
    mutate(lake_area_m2 = lake_area_km2*1000000) %>%
    select(-lake_area_km2) %>%
    na.omit(.) %>%
    mutate(ebullition_FO_mg_m2_d = ifelse(percent_ice_cover > 75, 0, mapply(FOA_model, mean_temp_k-273.15)),
           ebullition_Aben_mg_m2_d = ifelse(percent_ice_cover > 75, 0, mapply(Aben_model, mean_temp_k-273.15)),
           ebullition_LAC_FO = ebullition_FO_mg_m2_d*lake_area_m2*0.000001*0.6*30,
           ebullition_LAC_Aben = ebullition_Aben_mg_m2_d*lake_area_m2*0.000001*0.6*30) %>%
    write.table(., file = paste0("./output/global_ebullition_emission.csv"),
                   append = T,
                   row.names = F,
                   col.names = !file.exists("./output/global_ebullition_emission.csv"))
}

no_cores <- detectCores()-2
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
foreach(x=country) %dopar% emission_function(x)

#### Time check ####
e <- Sys.time()
t=e-s
print(t)


