# get the hybas ids and get the glcp data to link spatiallyand temporally. 

library(dplyr)
library(readr)
library(vroom)
library(doParallel)

country <- list.files(path = "./data/glcp_countries")
country <- gsub("\\..*", "", country)

analysis_function <- function(x){

obs <- read_csv("./data/organized_data_to_append/GLEE_data_with_hybas_ids.csv") %>%
  select(-published_data_product_DOI, -published_data_product, -primary_data_source_DOI) %>%
  filter(country == x) %>% 
  mutate(month = as.numeric(month))

obs_hybas_ids <- c(unique(obs$hybas_id))

glcp <- vroom::vroom(paste0("./data/glcp_countries/",x,".csv"), col_names = c("year", "month","country",
                                                                            "hybas_id","pop_sum","rel_hum",
                                                                            "sum_precip_mm", "mean_temp_k", 
                                                                            "mean_sw_wm2", "ice_cover_mean")) %>%
  filter(hybas_id %in% obs_hybas_ids) %>%
  select(-country, -ice_cover_mean) %>%
  group_by(hybas_id) %>%
  dplyr::mutate_at(vars(pop_sum),funs(imputeTS::na_interpolation(., option = "linear"))) %>%
  ungroup()

obs_new <- obs %>% left_join(., glcp, by = c("hybas_id", "year", "month")) %>%
  group_by(lat, lon, waterbody_type, year, month, country, continent, 
           primary_data_source) %>% 
  summarize_all(funs(mean), na.rm = T) %>%
  write.table(., file = paste0("./data/GLEE_data_with_GLCP_link.csv"),
              append = T,
              row.names = F,
              col.names = !file.exists("./data/GLEE_data_with_GLCP_link.csv"))
}

no_cores <- detectCores()-2
cl <- makeCluster(no_cores, type="FORK")
registerDoParallel(cl)
foreach(x=country) %dopar% analysis_function(x)

#### Time check ####
e <- Sys.time()
t=e-s
print(t)






