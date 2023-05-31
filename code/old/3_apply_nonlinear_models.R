library(minpack.lm)

ebu_base <- base %>% select(ch4_ebu, temp_for_model_K, surf_area_k, waterbody_id, month) %>%
  na.omit(.)%>%
  mutate(temp_for_model_C = temp_for_model_K-273.15)

# First-order Arhennius
monthly_ebu_FOA = nlsLM(ch4_ebu ~ A * exp(a * temp_for_model_C),
                     start = list(A = 0.12, a = 0.035),
                     data = ebu_base,
                     control = nls.lm.control(maxiter=1000))
summary(monthly_ebu_FOA) # get MSE value

# Modified First-order Arhennius (from Aben et al 2018 & Johnson et al 2022)
monthly_ebu_MFOA = nlsLM(ch4_ebu ~ A * a^(temp_for_model_C-20),
                     start = list(A = 100, a = 1.1),
                     data = ebu_base,
                     control = nls.lm.control(maxiter=1000))
summary(monthly_ebu_MFOA) # get MSE value

# Modified Second-Order Arhennius with surface area as limiting factor

monthly_ebu_MSOA = nlsLM(ch4_ebu ~ A * exp(a * temp_for_model_C) * exp(-(b * surf_area_k)),
                         start = list(A = 1.782e-07, a = 1, b=1),
                         data = ebu_base,
                         control = nls.lm.control(maxiter=1000))

summary(monthly_ebu_MSOA) # get MSE value
monthly_ebu_MFOA_function <- function(x,y) {7.473e-06 * exp(-7.962e-01*x) * exp(-(-1.408e-02*y))}


monthly_ebu_FOA_function <- function(x) {4.14072 * exp(0.15768*x)}
monthly_MFOA_function <- function(x) {100 * 1.1 ^ (x-20)}
jansen_ebu_FOA_function <- function(x)
aben_MFOA_function <- function(x) {100 * 1.1 ^ (x-20)}

ebu_function_output <- data.frame(x = -15:35,            # Create data for ggplot2
                                  values = c(
                                    monthly_ebu_FOA_function(-15:35),
                                    monthly_ebu_MFOA_function(-15:35),
                                    johnson_MFOA_function(-15:35)),
                                  model = rep(c(
                                    "First-Order Arhennius",
                                    "Modified First-Order Arhennius",
                                    "Aben's First-Order Arhennius"), each = 51))


ggplot(ebu_base, aes(x = temp_for_model_C, y = ch4_ebu))+
  geom_point()+
  geom_line(data = ebu_function_output, aes(x, values, group = model, color = model), lwd = 2)+
  scale_color_viridis(discrete = T, option = "C")+
  coord_cartesian(ylim=c(0, 3000))+
  theme_classic()



diff_base <- base %>% select(ch4_diff, temp_for_model_K) %>%
  na.omit(.)


monthly_diff_FKT <- nls(ch4_diff ~ A * exp(a * temp_for_model_K) 
                       , nls.control(maxiter=10000) , data = diff_base, start = list(A = 0.12, a = 0.035), trace = TRUE)

summary(monthly_diff_FKT) # get MSE value

temp_diff_function_mcclure <- function(x) {0.18630 * exp(0.01773*x)}

diff_function_output <- data.frame(x = 253:313,            # Create data for ggplot2
                                  values = c(
                                    temp_diff_function_mcclure(253:313)),
                                  model = rep(c(
                                    "New Calibration"), each = 61))

ggplot(diff_base, aes(x = temp_for_model_K, y = ch4_diff))+
  geom_point()+
  geom_line(data = diff_function_output, aes(x, values), lwd = 2)+
  theme_classic()

