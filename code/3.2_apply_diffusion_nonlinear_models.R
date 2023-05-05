library(minpack.lm)


### TEMPERATURE MODELS ###

### MONTHLY TIMESTEP ###
diff_base_temp <- base %>% select(ch4_diff, temp_for_model_K, waterbody_id, month) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15)

# First-order Arhennius
monthly_diff_FOA = nlsLM(ch4_ebu ~ A * exp(a * temp_for_model_C),
                        start = list(A = 0.023, a = 0.124),
                        data = diff_base_temp,
                        control = nls.lm.control(maxiter=1000))
summary(monthly_diff_FOA) # get MSE value

monthly_FOA_function <- function(x) {1.199e+02 * exp(4.895e-02*x)}
natchimuthu_FOA_function <- function(x) {0.023 * exp(0.124*x)}

diff_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
                                  values = c(
                                    monthly_FOA_function(-15:35),
                                    natchimuthu_FOA_function(-15:35)),
                                  model = rep(c(
                                    "First-Order (F-O) Arhennius",
                                    "Natchimuthu's F-O Arhennius"), each = 51))

ggplot(diff_base_temp, aes(x = temp_for_model_C, y = ch4_diff))+
  geom_point()+
  geom_line(data = diff_function_output, aes(x, values, group = model, color = model), lwd = 2)+
  scale_color_viridis(discrete = T, option = "C")+
  coord_cartesian(ylim=c(0, 2000))+
  theme_classic()



### YEARLY TIMESTEP ###

diff_base_yearly <- base %>% select(year, ch4_diff, waterbody_id, temp_for_model_K) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  group_by(waterbody_id, year)%>%
  summarize_all(funs(mean), na.rm = T)

# First-order Arhennius
yearly_diff_FOA = nlsLM(ch4_diff ~ A * exp(a * temp_for_model_C),
                       start = list(A = 0.12, a = 0.035),
                       data = diff_base_yearly,
                       control = nls.lm.control(maxiter=1000))
summary(yearly_diff_FOA) # get MSE value


yearly_FOA_function <- function(x) {18.05549 * exp(0.03913*x)}
natchimuthu_FOA_function <- function(x) {0.023 * exp(0.124*x)}

diff_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
                                   values = c(
                                     yearly_FOA_function(-15:35),
                                     natchimuthu_FOA_function(-15:35)),
                                   model = rep(c(
                                     "First-Order (F-O) Arhennius",
                                     "Natchimuthu's F-O Arhennius"), each = 51))

ggplot(diff_base_yearly, aes(x = temp_for_model_C, y = ch4_diff))+
  geom_point()+
  geom_line(data = diff_function_output, aes(x, values, group = model, color = model), lwd = 2)+
  scale_color_viridis(discrete = T, option = "C")+
  coord_cartesian(ylim=c(0, 3000))+
  theme_classic()

### SITE MEAN ###

diff_base_site <- base %>% select(ch4_diff, waterbody_id, temp_for_model_K) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15) %>%
  group_by(waterbody_id)%>%
  summarize_all(funs(mean), na.rm = T)

# First-order Arhennius
site_diff_FOA = nlsLM(ch4_diff ~ A * exp(a * temp_for_model_C),
                     start = list(A = 0.12, a = 0.035),
                     data = diff_base_site,
                     control = nls.lm.control(maxiter=1000))
summary(site_diff_FOA) # get MSE value

site_FOA_function <- function(x) {19.37369 * exp(0.03704*x)}
natchimuthu_FOA_function <- function(x) {0.023 * exp(0.124*x)}

diff_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
                                   values = c(
                                     yearly_FOA_function(-15:35),
                                     natchimuthu_FOA_function(-15:35)),
                                   model = rep(c(
                                     "First-Order (F-O) Arhennius",
                                     "Natchimuthu's F-O Arhennius"), each = 51))

ggplot(diff_base_site, aes(x = temp_for_model_C, y = ch4_diff))+
  geom_point()+
  geom_line(data = diff_function_output, aes(x, values, group = model, color = model), lwd = 2)+
  scale_color_viridis(discrete = T, option = "C")+
  coord_cartesian(ylim=c(0, 3000))+
  theme_classic()


### LAKE AREA MODELS ###

### LAKE AREA ALL DATA ###

diff_base_area <- base %>% select(ch4_diff, waterbody_id, month, surf_area_k) %>%
  na.omit(.) %>%
  mutate(ch4_diff_scaled = scale(ch4_diff),
         surf_area_k_scaled = scale(surf_area_k))

area_monthly_model <- lm(ch4_diff_scaled ~ surf_area_k_scaled, data = diff_base_area)
summary(area_monthly_model)




