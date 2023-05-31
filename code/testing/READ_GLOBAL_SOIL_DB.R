library(ncdf4)
library(sf)

# set path and filename
ncpath <- "/Users/ryanmcclure/Documents/GLEE2/data/HWSD_1247/data/"
ncname <- "AWT_T_SOC"  
ncfname <- paste(ncpath, ncname, ".nc4", sep="")

ncin <- nc_open(ncfname)
print(ncin)

lon <- ncvar_get(ncin,"lon")
lat <- ncvar_get(ncin,"lat")
T_OC <- ncvar_get(ncin,"SUM_t_c_12")


# get temperature
T_OC_array <- ncvar_get(ncin,"SUM_t_c_12")
dlname <- ncatt_get(ncin,"SUM_t_c_12","long_name")
dunits <- ncatt_get(ncin,"SUM_t_c_12","units")
fillvalue <- ncatt_get(ncin,"SUM_t_c_12","_FillValue")
dim(T_OC_array)

title <- ncatt_get(ncin,0,"title")
institution <- ncatt_get(ncin,0,"institution")
datasource <- ncatt_get(ncin,0,"source")
references <- ncatt_get(ncin,0,"references")
history <- ncatt_get(ncin,0,"history")
Conventions <- ncatt_get(ncin,0,"Conventions")

library(lattice)
library(RColorBrewer)

lonlat <- as.matrix(expand.grid(lon,lat))
T_OC_vec <- as.vector(T_OC_array)

TOC_df <- data.frame(cbind(lonlat,T_OC_vec))
names(TOC_df) <- c("lon","lat",paste("T_OC", sep="_"))
head(na.omit(TOC_df), 10)

TOC_df <- TOC_df %>%
st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext")

base2 <- base %>%  
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform("+proj=eqearth +wktext") %>%
  st_join(., TOC_df, join = st_nearest_feature)


base2 <- read_csv("./data/GLEE_data_with_GLCP_HWSD_link.csv") %>%
  mutate(tot_sampling_events = num_months_sampled * num_sites_sampled) %>%
  group_by(geometry) %>%
  mutate(waterbody_id = cur_group_id()) %>%
  ungroup() %>%
  mutate(temp_for_model_K = ifelse(is.na(mean_temp_k), effective_obs_wtemp_k, mean_temp_k))

ebu_base_temp_area <- base2 %>% select(ch4_ebu, temp_for_model_K, waterbody_id, month, surf_area_k) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15)

monthly_ebu_SOA_area_threshold = nlsLM(ch4_ebu ~ A * (exp(a * temp_for_model_C) * (surf_area_k/k+surf_area_k)),
                                       start = list(A = 0.12, a = 0.035, k = 100),
                                       data = ebu_base_temp_area,
                                       control = nls.lm.control(maxiter=1000))
summary(monthly_ebu_SOA_area_threshold) # get MSE value

dat <- as.data.frame(cbind(ebu_base_temp_area$ch4_ebu, predict(monthly_ebu_SOA_area_threshold)))
second_order_area_limiter_NSE <- NSE(dat$V2, dat$V1)


ebu_base_temp_soil <- base2 %>% select(ch4_ebu, temp_for_model_K, waterbody_id, month, T_OC) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15)

monthly_ebu_SOA_area_threshold = nlsLM(ch4_ebu ~ A * exp(a * temp_for_model_C) * (T_OC/k+T_OC),
                                       start = list(A = 0.12, a = 0.035, k = 100),
                                       data = ebu_base_temp_soil,
                                       control = nls.lm.control(maxiter=1000))
summary(monthly_ebu_SOA_area_threshold) # get MSE value

dat <- as.data.frame(cbind(ebu_base_temp_soil$ch4_ebu, predict(monthly_ebu_SOA_area_threshold)))
second_order_area_limiter_NSE <- NSE(dat$V2, dat$V1)
