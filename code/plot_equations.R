


p1 <- ggplot(ebu_base_temp_lake, aes(x = temp_for_model_C, y = ch4_ebu))+
  geom_point(size = 4, pch = 21, fill = "purple")+
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Ebullition Rate")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))
ggsave("./output/lake1.png", width = 12, height = 8)


p2 <- ggplot(baseline_lake_ebu, aes(x = temp_for_model_C, y = fit))+
  geom_point(data = ebu_base_temp_lake, aes(x = temp_for_model_C, y = ch4_ebu), size = 4, pch = 21, fill = "purple")+
  geom_line(data = baseline_lake_ebu, aes(x = temp_for_model_C, y = fit), color = "green4", lwd = 2) +
  geom_ribbon(data = baseline_lake_ebu, aes(ymin = low_CI, ymax = high_CI), fill = "green4", color = NA,alpha = 0.2)+
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Ebullition Rate")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))
ggsave("./output/lake_base.png", width = 12, height = 8)

p3 <- ggplot(time_lake_ebu, aes(x = temp_for_model_C, y = fit))+
  geom_point(data = ebu_base_temp_lake, aes(x = temp_for_model_C, y = ch4_ebu), size = 4, pch = 21, fill = "purple")+
  geom_point(data = ebu_base_temp_lake, aes(x = temp_for_model_C, y = ebu_time_low), size = 2, pch = 21, fill = "grey70")+
  geom_point(data = ebu_base_temp_lake, aes(x = temp_for_model_C, y = ebu_time_high), size = 2, pch = 21, fill = "grey1")+
  geom_line(data = time_lake_ebu, aes(x = temp_for_model_C, y = fit, group = scenario), color = "blue", lwd = 2) +
  geom_ribbon(data = time_lake_ebu, aes(ymin = low_CI, ymax = high_CI, group = scenario), fill = "blue", color = NA,alpha = 0.2)+
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Ebullition (mg CH4 m-2 d-1)")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))
ggsave("./output/lake_time_error.png", width = 12, height = 8)

p4 <- ggplot(space_lake_ebu, aes(x = temp_for_model_C, y = fit))+
  geom_point(data = ebu_base_temp_lake, aes(x = temp_for_model_C, y = ch4_ebu), size = 4, pch = 21, fill = "purple")+
  geom_point(data = ebu_base_temp_lake, aes(x = temp_for_model_C, y = ebu_space_low), size = 2, pch = 21, fill = "grey70")+
  geom_point(data = ebu_base_temp_lake, aes(x = temp_for_model_C, y = ebu_space_high), size = 2, pch = 21, fill = "grey1")+
  geom_line(data = space_lake_ebu, aes(x = temp_for_model_C, y = fit, group = scenario), color = "red", lwd = 2) +
  geom_ribbon(data = space_lake_ebu, aes(ymin = low_CI, ymax = high_CI, group = scenario), fill = "red", color = NA,alpha = 0.2)+
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Ebullition (mg CH4 m-2 d-1)")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))
ggsave("./output/lake_space_error.png", width = 12, height = 8)



# p5 <- ggplot(error_ebu_lake, aes(x = temp_for_model_C, y = ch4_ebu))+
#   geom_point(size = 4, pch = 21, fill = "purple")+
#   geom_line(data = ebullition_outputs_lake_param, aes(x, values, group = model),color = "orange", lwd = 2)+
#   scale_color_viridis(discrete = T, option = "C")+
#   coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
#   ylab("Ebullition Emission Rate")+
#   xlab("Temperature (C)")+
#   theme_classic()+
#   theme(axis.title = element_text(size = 15, color = "black"),
#         axis.text = element_text(size = 15, color = "black"))
# ggsave("./output/lake_param_error.png", width = 12, height = 8)
# 
# p6 <- ggplot(error_ebu_lake, aes(x = temp_for_model_C, y = ch4_ebu))+
#   geom_point(size = 4, pch = 21, fill = "purple")+
#   geom_line(data = ebullition_outputs_lake_model, aes(x, values, group = model),color = "grey50", lwd = 2)+
#   scale_color_viridis(discrete = T, option = "C")+
#   coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
#   ylab("Ebullition Emission Rate")+
#   xlab("Temperature (C)")+
#   theme_classic()+
#   theme(axis.title = element_text(size = 15, color = "black"),
#         axis.text = element_text(size = 15, color = "black"))
# ggsave("./output/lake_model_error.png", width = 12, height = 8)
# 
# p7 <- ggplot(error_ebu_lake, aes(x = temp_for_model_C, y = ch4_ebu))+
#   geom_point(size = 4, pch = 21, fill = "purple")+
#   coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
#   ylab("Ebullition Emission Rate")+
#   xlab("Temperature (C)")+
#   theme_classic()+
#   theme(axis.title = element_text(size = 15, color = "black"),
#         axis.text = element_text(size = 15, color = "black"))
# ggsave("./output/lake2.png", width = 12, height = 8)
# 
# p8 <- ggplot(error_ebu_lake, aes(x = temp_for_model_C, y = ch4_ebu))+
#   geom_point(size = 4, pch = 21, fill = "purple")+
#   geom_point(data = error_ebu_res, aes(temp_for_model_C, ch4_ebu), size = 4, pch = 21, fill = "cyan")+
#   scale_color_viridis(discrete = T, option = "C")+
#   coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
#   ylab("Ebullition Emission Rate")+
#   xlab("Temperature (C)")+
#   theme_classic()+
#   theme(axis.title = element_text(size = 15, color = "black"),
#         axis.text = element_text(size = 15, color = "black"))
# ggsave("./output/lake_res.png", width = 12, height = 8)


p9 <- ggplot(ebu_base_temp_res, aes(x = temp_for_model_C, y = ch4_ebu))+
  geom_point(size = 4, pch = 21, fill = "cyan")+
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Ebullition Rate")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))
ggsave("./output/res.png", width = 12, height = 8)

p10 <- ggplot(baseline_res_ebu, aes(x = temp_for_model_C, y = fit))+
  geom_point(data = ebu_base_temp_res, aes(x = temp_for_model_C, y = ch4_ebu), size = 4, pch = 21, fill = "cyan")+
  geom_line(data = baseline_res_ebu, aes(x = temp_for_model_C, y = fit), color = "green4", lwd = 2) +
  geom_ribbon(data = baseline_res_ebu, aes(ymin = low_CI, ymax = high_CI), fill = "green4", color = NA,alpha = 0.2)+
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Ebullition Rate")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))
ggsave("./output/res_base.png", width = 12, height = 8)

p11 <- ggplot(time_res_ebu, aes(x = temp_for_model_C, y = fit))+
  geom_point(data = ebu_base_temp_res, aes(x = temp_for_model_C, y = ch4_ebu), size = 4, pch = 21, fill = "cyan")+
  geom_point(data = ebu_base_temp_res, aes(x = temp_for_model_C, y = ebu_time_low), size = 2, pch = 21, fill = "grey70")+
  geom_point(data = ebu_base_temp_res, aes(x = temp_for_model_C, y = ebu_time_high), size = 2, pch = 21, fill = "grey1")+
  geom_line(data = time_res_ebu, aes(x = temp_for_model_C, y = fit, group = scenario), color = "blue", lwd = 2) +
  geom_ribbon(data = time_res_ebu, aes(ymin = low_CI, ymax = high_CI, group = scenario), fill = "blue", color = NA,alpha = 0.1)+
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Ebullition (mg CH4 m-2 d-1)")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))
ggsave("./output/res_time_error.png", width = 12, height = 8)

p12 <-ggplot(space_res_ebu, aes(x = temp_for_model_C, y = fit))+
  geom_point(data = ebu_base_temp_res, aes(x = temp_for_model_C, y = ch4_ebu), size = 4, pch = 21, fill = "cyan")+
  geom_point(data = ebu_base_temp_res, aes(x = temp_for_model_C, y = ebu_space_low), size = 2, pch = 21, fill = "grey70")+
  geom_point(data = ebu_base_temp_res, aes(x = temp_for_model_C, y = ebu_space_high), size = 2, pch = 21, fill = "grey1")+
  geom_line(data = space_res_ebu, aes(x = temp_for_model_C, y = fit, group = scenario), color = "red", lwd = 2) +
  geom_ribbon(data = space_res_ebu, aes(ymin = low_CI, ymax = high_CI, group = scenario), fill = "red", color = NA,alpha = 0.2)+
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Ebullition (mg CH4 m-2 d-1)")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))
ggsave("./output/res_space_error.png", width = 12, height = 8)

p13 <- ggplot(error_ebu_res, aes(x = temp_for_model_C, y = ch4_ebu))+
  geom_point(size = 3, pch = 21, fill = "cyan")+
  geom_line(data = ebullition_outputs_res_param, aes(x, values, group = model), color = "orange", lwd = 2)+
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Ebullition Emission Rate")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))
ggsave("./output/res_param_error.png", width = 12, height = 8)

p14 <- ggplot(error_ebu_res, aes(x = temp_for_model_C, y = ch4_ebu))+
  geom_point(size = 3, pch = 21, fill = "cyan")+
  geom_line(data = ebullition_outputs_res_model, aes(x, values, group = model), color = "grey50", lwd = 2)+
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Ebullition Emission Rate")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))

ggsave("./output/res_model_error.png", width = 12, height = 8)



p15 <- ggplot(error_ebu_lake, aes(x = temp_for_model_C, y = ch4_ebu))+
  geom_point(size = 3, pch = 21, fill = "purple")+
  geom_line(data = ebullition_outputs_lake_base, aes(x, values, group = model), color = "green3", lwd = 1) +
  geom_line(data = ebullition_outputs_lake_time, aes(x, values, group = model),color = "blue", lwd = 1) +
  geom_line(data = ebullition_outputs_lake_space, aes(x, values, group = model),color = "red", lwd = 1) +
  geom_line(data = ebullition_outputs_lake_param, aes(x, values, group = model),color = "orange", lwd = 1)+
  geom_line(data = ebullition_outputs_lake_model, aes(x, values, group = model), color = "grey50", lwd = 1) +
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Ebullition Emission Rate")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))

ggsave("./output/lake_all.png", width = 12, height = 6)

p16 <- ggplot(error_ebu_res, aes(x = temp_for_model_C, y = ch4_ebu))+
  geom_point(size = 3, pch = 21, fill = "cyan")+
  geom_line(data = ebullition_outputs_res_base, aes(x, values, group = model), color = "green3", lwd = 1) +
  geom_line(data = ebullition_outputs_res_time, aes(x, values, group = model),color = "blue", lwd = 1) +
  geom_line(data = ebullition_outputs_res_space, aes(x, values, group = model),color = "red", lwd = 1) +
  geom_line(data = ebullition_outputs_res_param, aes(x, values, group = model),color = "orange", lwd = 1)+
  geom_line(data = ebullition_outputs_res_model, aes(x, values, group = model), color = "grey50", lwd = 1) +
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Ebullition Emission Rate")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))

ggsave("./output/res_all.png", width = 12, height = 6)







source("./code/diffusion_functions.R")

diffusion_outputs <- data.frame(x = -5:35,           
                                 values = c(
                                   FOA_mean_lake(-5:35),
                                   FOA_time_low_lake(-5:35),
                                   FOA_time_high_lake(-5:35),
                                   FOA_space_low_lake(-5:35),
                                   FOA_space_high_lake(-5:35),
                                   FOA_coefficient_low_lake(-5:35),
                                   FOA_coefficient_high_lake(-5:35),
                                   FOA_mean_error_low_lake(-5:35),
                                   FOA_mean_error_high_lake(-5:35),
                                   FOA_mean_res(-5:35),
                                   FOA_time_low_res(-5:35),
                                   FOA_time_high_res(-5:35),
                                   FOA_space_low_res(-5:35),
                                   FOA_space_high_res(-5:35),
                                   FOA_coefficient_low_res(-5:35),
                                   FOA_coefficient_high_res(-5:35),
                                   FOA_mean_error_low_res(-5:35),
                                   FOA_mean_error_high_res(-5:35)),
                                 model = rep(c(
                                   "Baseline (lake)",
                                   "Temporal -1 S.D. (lake)",
                                   "Temporal +1 S.D. (lake)",
                                   "Spatial -1 S.D. (lake)",
                                   "Spatial +1 S.D. (lake)",
                                   "Coefficients -1 S.D. (lake)",
                                   "Coefficients +1 S.D. (lake)",
                                   "Model -1 S.D. (lake)",
                                   "Model +1 S.D. (lake)",
                                   "Baseline (res)",
                                   "Temporal -1 S.D. (res)",
                                   "Temporal +1 S.D. (res)",
                                   "Spatial -1 S.D. (res)",
                                   "Spatial +1 S.D. (res)",
                                   "Coefficients -1 S.D. (res)",
                                   "Coefficients +1 S.D. (res)",
                                   "Model -1 S.D. (res)",
                                   "Model +1 S.D. (res)"), 
                                   each = 41))

error_diff_lake <- error_diff %>% filter(waterbody_type == "lake")
error_diff_res <- error_diff %>% filter(waterbody_type == "reservoir")


diffusion_outputs_lake <- diffusion_outputs %>% filter(grepl('(lake)', model))
diffusion_outputs_lake_base <- diffusion_outputs_lake %>% filter(grepl('Baseline', model))
diffusion_outputs_lake_time <- diffusion_outputs_lake %>% filter(grepl('Temporal', model))
diffusion_outputs_lake_space <- diffusion_outputs_lake %>% filter(grepl('Spatial', model))
diffusion_outputs_lake_param <- diffusion_outputs_lake %>% filter(grepl('Coefficient', model))
diffusion_outputs_lake_model <- diffusion_outputs_lake %>% filter(grepl('Model', model)) %>%
  mutate(values = ifelse(values<0,0,values))


diffusion_outputs_res <- diffusion_outputs %>% filter(grepl('(res)', model))

diffusion_outputs_res_base <- diffusion_outputs_res %>% filter(grepl('Baseline', model))
diffusion_outputs_res_time <- diffusion_outputs_res %>% filter(grepl('Temporal', model))
diffusion_outputs_res_space <- diffusion_outputs_res %>% filter(grepl('Spatial', model))
diffusion_outputs_res_param <- diffusion_outputs_res %>% filter(grepl('Coefficient', model))
diffusion_outputs_res_model <- diffusion_outputs_res %>% filter(grepl('Model', model))%>%
  mutate(values = ifelse(values<0,0,values))

p15 <- ggplot(error_diff_lake, aes(x = temp_for_model_C, y = ch4_diff))+
  geom_point(size = 3, pch = 21, fill = "purple")+
  geom_line(data = diffusion_outputs_lake_base, aes(x, values, group = model), color = "green3", lwd = 1) +
  geom_line(data = diffusion_outputs_lake_time, aes(x, values, group = model),color = "blue", lwd = 1) +
  geom_line(data = diffusion_outputs_lake_space, aes(x, values, group = model),color = "red", lwd = 1) +
  geom_line(data = diffusion_outputs_lake_param, aes(x, values, group = model),color = "orange", lwd = 1)+
  geom_line(data = diffusion_outputs_lake_model, aes(x, values, group = model), color = "grey50", lwd = 1) +
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Diffusion Emission Rate")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))

ggsave("./output/diff_lake_all.png", width = 12, height = 6)

p16 <- ggplot(error_diff_res, aes(x = temp_for_model_C, y = ch4_diff))+
  geom_point(size = 3, pch = 21, fill = "cyan")+
  geom_line(data = diffusion_outputs_res_base, aes(x, values, group = model), color = "green3", lwd = 1) +
  geom_line(data = diffusion_outputs_res_time, aes(x, values, group = model),color = "blue", lwd = 1) +
  geom_line(data = diffusion_outputs_res_space, aes(x, values, group = model),color = "red", lwd = 1) +
  geom_line(data = diffusion_outputs_res_param, aes(x, values, group = model),color = "orange", lwd = 1)+
  geom_line(data = diffusion_outputs_res_model, aes(x, values, group = model), color = "grey50", lwd = 1) +
  coord_cartesian(ylim=c(0, 1500),xlim=c(0, 35))+
  ylab("Diffusion Emission Rate")+
  xlab("Temperature (C)")+
  theme_classic()+
  theme(axis.title = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 15, color = "black"))

ggsave("./output/diff_res_all.png", width = 12, height = 6)
