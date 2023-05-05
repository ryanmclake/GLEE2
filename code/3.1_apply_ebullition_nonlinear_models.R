library(minpack.lm)
library(hydroGOF)

### TEMPERATURE MODELS ###

### MONTHLY TIMESTEP ###
ebu_base_temp <- base %>% select(ch4_ebu, temp_for_model_K, waterbody_id, month) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15)

# First-order Arhennius
monthly_ebu_SOA_temp_threshold = nlsLM(ch4_ebu ~ A * exp(a * temp_for_model_C - b * temp_for_model_C^2),
                        start = list(A = 0.12, a = 0.035, b = 0.0001),
                        data = ebu_base_temp,
                        control = nls.lm.control(maxiter=1000))
summary(monthly_ebu_SOA_temp_threshold) # get MSE value

dat <- as.data.frame(cbind(ebu_base_temp$ch4_ebu, predict(monthly_ebu_SOA_temp_threshold)))

second_order_temp_limiter_NSE <- NSE(dat$V2, dat$V1)

ebu_base_temp_area <- base %>% select(ch4_ebu, temp_for_model_K, waterbody_id, month, surf_area_k) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15)

monthly_ebu_SOA_area_threshold = nlsLM(ch4_ebu ~ A * exp(a * temp_for_model_C - b * surf_area_k^2),
                        start = list(A = 0.12, a = 0.035, b = 0.0001),
                        data = ebu_base_temp_area,
                        control = nls.lm.control(maxiter=1000))
summary(monthly_ebu_SOA_area_threshold) # get MSE value

dat <- as.data.frame(cbind(ebu_base_temp_area$ch4_ebu, predict(monthly_ebu_SOA_area_threshold)))
second_order_area_limiter_NSE <- NSE(dat$V2, dat$V1)

# Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
monthly_ebu_MFOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
                         start = list(A = 100, a = 1.1),
                         data = ebu_base_temp,
                         control = nls.lm.control(maxiter=1000))
summary(monthly_ebu_MFOA) # get MSE value

dat <- as.data.frame(cbind(ebu_base_temp$ch4_ebu, predict(monthly_ebu_MFOA)))

first_order_temp_rmse <- NSE(dat$V2, dat$V1)

# Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
monthly_ebu_MFOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
                         start = list(A = 100, a = 1.1),
                         data = ebu_base_temp,
                         control = nls.lm.control(maxiter=1000))
summary(monthly_ebu_MFOA) # get MSE value



monthly_ebu_SOA_function <- function(x) {30.136857 * exp(0.016959*x - (-0.003697*x^2))}
monthly_MFOA_function <- function(x) {96.98125 * 1.17077 ^ (x-20)}
aben_MFOA_function <- function(x) {100 * 1.1 ^ (x-20)}


ebu_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
                                  values = c(
                                   monthly_MFOA_function(-15:35),
                                   aben_MFOA_function(-15:35),
                                   monthly_ebu_SOA_function(-15:35)),
                                  model = rep(c(
                                    "First-Order Arhennius",
                                    "Aben's Arhennius",
                                    "Second-Order Exponential"), each = 51))

ggplot(ebu_base_temp, aes(x = temp_for_model_C, y = ch4_ebu))+
  geom_point()+
  geom_line(data = ebu_function_output, aes(x, values, group = model, color = model), lwd = 2)+
  scale_color_viridis(discrete = T, option = "C")+
  coord_cartesian(ylim=c(0, 3000))+
  theme_classic()



### YEARLY TIMESTEP ###

ebu_base_yearly <- base %>% select(year, ch4_ebu, waterbody_id, temp_for_model_K) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  group_by(waterbody_id, year)%>%
  summarize_all(funs(mean), na.rm = T)

# First-order Arhennius
yearly_ebu_FOA = nlsLM(ch4_ebu ~ A * exp(a * temp_for_model_C),
                        start = list(A = 0.12, a = 0.035),
                        data = ebu_base_yearly,
                        control = nls.lm.control(maxiter=1000))
summary(yearly_ebu_FOA) # get MSE value

# Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
yearly_ebu_MFOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
                         start = list(A = 100, a = 1.1),
                         data = ebu_base_yearly,
                         control = nls.lm.control(maxiter=1000))
summary(yearly_ebu_MFOA) # get MSE value

# Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
yearly_ebu_MFOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
                         start = list(A = 100, a = 1.1),
                         data = ebu_base_yearly,
                         control = nls.lm.control(maxiter=1000))
summary(yearly_ebu_MFOA) # get MSE value

yearly_FOA_function <- function(x) {13.71050 * exp(0.08813*x)}
yearly_MFOA_function <- function(x) {79.89646 * 1.09213 ^ (x-20)}
aben_MFOA_function <- function(x) {100 * 1.1 ^ (x-20)}

ebu_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
                                  values = c(
                                    yearly_FOA_function(-15:35),
                                    yearly_MFOA_function(-15:35),
                                    aben_MFOA_function(-15:35)),
                                  model = rep(c(
                                    "First-Order (F-O) Arhennius",
                                    "Modified F-O Arhennius",
                                    "Aben's Modified F-O Arhennius"), each = 51))

ggplot(ebu_base_yearly, aes(x = temp_for_model_C, y = ch4_ebu))+
  geom_point()+
  geom_line(data = ebu_function_output, aes(x, values, group = model, color = model), lwd = 2)+
  scale_color_viridis(discrete = T, option = "C")+
  coord_cartesian(ylim=c(0, 3000))+
  theme_classic()

### SITE MEAN ###

ebu_base_site <- base %>% select(ch4_ebu, waterbody_id, temp_for_model_K) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  group_by(waterbody_id)%>%
  summarize_all(funs(mean), na.rm = T)

# First-order Arhennius
site_ebu_FOA = nlsLM(ch4_ebu ~ A * exp(a * temp_for_model_C),
                       start = list(A = 0.12, a = 0.035),
                       data = ebu_base_site,
                       control = nls.lm.control(maxiter=1000))
summary(site_ebu_FOA) # get MSE value

# Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
site_ebu_MFOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
                        start = list(A = 100, a = 1.1),
                        data = ebu_base_site,
                        control = nls.lm.control(maxiter=1000))
summary(site_ebu_MFOA) # get MSE value

site_FOA_function <- function(x) {19.18062 * exp(0.06981*x)}
site_MFOA_function <- function(x) {77.48887 * 1.07230 ^ (x-20)}
aben_MFOA_function <- function(x) {100 * 1.1 ^ (x-20)}

ebu_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
                                  values = c(
                                    site_FOA_function(-15:35),
                                    site_MFOA_function(-15:35),
                                    aben_MFOA_function(-15:35)),
                                  model = rep(c(
                                    "First-Order (F-O) Arhennius",
                                    "Modified F-O Arhennius",
                                    "Aben's Modified F-O Arhennius"), each = 51))

ggplot(ebu_base_site, aes(x = temp_for_model_C, y = ch4_ebu))+
  geom_point()+
  geom_line(data = ebu_function_output, aes(x, values, group = model, color = model), lwd = 2)+
  scale_color_viridis(discrete = T, option = "C")+
  coord_cartesian(ylim=c(0, 3000))+
  theme_classic()


### LAKE AREA MODELS ###

### LAKE AREA ALL DATA ###

ebu_base_area <- base %>% select(ch4_ebu, waterbody_id, month, surf_area_k) %>%
  na.omit(.) %>%
  mutate(ch4_ebu_scaled = scale(ch4_ebu),
         surf_area_k_scaled = scale(surf_area_k))

area_monthly_model <- lm(ch4_ebu_scaled ~ surf_area_k_scaled, data = ebu_base_area)
summary(area_monthly_model)




ebu_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
                                  values = c(
                                    monthly_MFOA_function(-15:35),
                                    yearly_MFOA_function(-15:35),
                                    site_MFOA_function(-15:35),
                                    aben_MFOA_function(-15:35)),
                                  model = rep(c(
                                    "Monthly First-Order Arhennius",
                                    "Annual First-Order Arhennius",
                                    "Site-Level First-Order Arhennius",
                                    "Aben's Modified F-O Arhennius"), each = 51))

ggplot(ebu_base_temp, aes(x = temp_for_model_C, y = ch4_ebu))+
  geom_point()+
  geom_line(data = ebu_function_output, aes(x, values, group = model, color = model), lwd = 2)+
  scale_color_viridis(discrete = T, option = "C")+
  coord_cartesian(ylim=c(0, 3000))+
  theme_classic()






### Lake and reservoir partition ###

ebu_base <- base %>% select(ch4_ebu, temp_for_model_K, waterbody_type, month) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15,
         waterbody_type = ifelse(waterbody_type == "pond", "lake", waterbody_type)) 

ebu_base_reservoir <- base %>% select(ch4_ebu, temp_for_model_K, waterbody_type, month) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  filter(waterbody_type == "reservoir")

ebu_base_lake <- base %>% select(ch4_ebu, temp_for_model_K, waterbody_type, month) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  filter(waterbody_type == "lake")

# Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
reservoir_ebu_FOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
                         start = list(A = 100, a = 1.1),
                         data = ebu_base_reservoir,
                         control = nls.lm.control(maxiter=1000))
summary(reservoir_ebu_FOA) # get MSE value



# Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
lake_ebu_FOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
                         start = list(A = 100, a = 1.1),
                         data = ebu_base_lake,
                         control = nls.lm.control(maxiter=1000))
summary(lake_ebu_FOA) # get MSE value


lake_FOA_function <- function(x) {77.07198 * 1.04613 ^ (x-20)}
reservoir_MFOA_function <- function(x) {10.3765 * 1.6844 ^ (x-20)}
aben_MFOA_function <- function(x) {100 * 1.1 ^ (x-20)}

ebu_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
                                  values = c(
                                    lake_FOA_function(-15:35),
                                    reservoir_MFOA_function(-15:35),
                                    aben_MFOA_function(-15:35)),
                                  model = rep(c(
                                    "First-Order Arhennius (Lake)",
                                    "First-Order Arhennius (Reservoir)",
                                    "Aben's Arhennius"), each = 51))

ggplot(ebu_base, aes(x = temp_for_model_C, y = ch4_ebu))+
  geom_point(aes(shape = waterbody_type), size = 3)+
  geom_line(data = ebu_function_output, aes(x, values, group = model, color = model), lwd = 2)+
  scale_color_viridis(discrete = T, option = "C")+
  coord_cartesian(ylim=c(0, 3000))+
  theme_classic()



function1 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 0.001^2))} 
function2 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 0.01^2))} 
function3 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 0.1^2))} 
function4 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 1^2))} 
function5 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 10^2))} 
function6 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 100^2))} 
function7 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 1000^2))} 
function8 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 2000^2))} 
function9 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 3000^2))} 
function10 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 4000^2))} 
function11 <- function(x) {(2.385e-06 * exp(6.975e-01 * x - 3.765e-07 * 5000^2))} 




function_output <- data.frame(x = -15:35,            # Create data for ggplot2
                              values = c(function1(-15:35),
                                         function2(-15:35),
                                         function3(-15:35),
                                         function4(-15:35),
                                         function5(-15:35),
                                         function6(-15:35),
                                         function7(-15:35),
                                         function8(-15:35),
                                         function9(-15:35),
                                         function10(-15:35),
                                         function11(-15:35)),
                              surface_area = rep(c(0.001,
                                                0.01,
                                                0.1,
                                                1,
                                                10,
                                                100,
                                                1000,
                                                2000,
                                                3000,
                                                4000,
                                                5000), each = 51)) %>%
  rename(`Lake Surface Area (km2)` = surface_area)



function_plot <- ggplot(function_output,
                        aes(x, values, group = `Lake Surface Area (km2)`, col = `Lake Surface Area (km2)`)) +
  geom_line()+
  scale_color_viridis(discrete = F, option = "C")+
  geom_point(data = ebu_base_temp_area, aes(x = temp_for_model_C, y = ch4_ebu), inherit.aes = F)+
  theme_light()+
  ylim(c(0,3000))+
  labs(x = "Water Temperature (C)", y = "Ebullition Rate")

ggsave(function_plot, path = ".",
       filename = "./figures/temp_do_elevation_concept_plot.jpg",
       width = 8, height = 6, device='jpg', dpi=1000)
