
library(tidyverse)

dat <- read_csv("./global_predictions.csv") %>%
  rename(`Variability Source` = group)


summary(dat)

cols <- c("a: Past" = "black",
          "b: Baseline" = "green", 
          "c: Time-low" = "lightblue1", 
          "d: Time-high" = "blue2",
          "e: Space-low" = "pink2", 
          "f: Space-high" = "red2",
          "g: Parameter-low" = "darkgoldenrod1", 
          "h: Parameter-high" = "darkorange3",
          "i: Model-low" = "grey70", 
          "j: Model-high" = "grey40")



p <- ggplot(dat, aes(x = area_weighted_flux, y = author, fill = `Variability Source`))+

ggplot(dat, aes(x = area_weighted_flux, y = author, fill = group))+
  scale_y_discrete(limits=rev)+
  theme_classic() + 
  scale_colour_manual(values = cols)+
  scale_fill_manual(values = cols)+
  geom_segment(aes(x = 72.85 , y = 10, xend = 72.85, yend = 5.3), lty = "dotted", lwd = 0.4, color = "black")+
  geom_segment(aes(x = 33.61 , y = 5, xend = 33.61, yend = 0), lty = "dashed", lwd = 0.4, color = "grey60")+
  geom_segment(aes(x = 33.61 , y = 5.2, xend = 72.85, yend = 5.2), lty = "solid", lwd = 0.2, color = "black", 
               arrow = grid::arrow(length = unit(0.1,"cm")))+
  annotate("text", x=50, y=5.3, label = "39", size=2)+
  geom_segment(aes(x = 22.94076 , y = 4.2, xend = 30, yend = 4.2), lty = "solid", lwd = 0.2, color = "black", #time
               arrow = grid::arrow(length = unit(0.1,"cm")))+
  geom_segment(aes(x = 162.6627, y = 4.2, xend = 36, yend = 4.2), lty = "solid", lwd = 0.2, color = "black", #time
               arrow = grid::arrow(length = unit(0.1,"cm")))+
  annotate("text", x=15, y=4.3, label = "11", size=2)+
  annotate("text", x=75, y=4.3, label = "130", size=2)+
  geom_segment(aes(x = 1.899455 , y = 3.2, xend = 30, yend = 3.2), lty = "solid", lwd = 0.2, color = "black", #space
               arrow = grid::arrow(length = unit(0.1,"cm")))+
  geom_segment(aes(x = 91.0399, y = 3.2, xend = 36, yend = 3.2), lty = "solid", lwd = 0.2, color = "black", #space
               arrow = grid::arrow(length = unit(0.1,"cm")))+
  annotate("text", x=10, y=3.3, label = "32", size=2)+
  annotate("text", x=65, y=3.3, label = "57", size=2)+
  geom_segment(aes(x =  31.797, y = 2.2, xend = 30, yend = 2.2), lty = "solid", lwd = 0.2, color = "black", #parameter
               arrow = grid::arrow(length = unit(0.1,"cm")))+
  geom_segment(aes(x = 216.0961, y = 2.2, xend = 36, yend = 2.2), lty = "solid", lwd = 0.2, color = "black", #parameter
               arrow = grid::arrow(length = unit(0.1,"cm")))+
  annotate("text", x=28, y=2.3, label = "2", size=2)+
  annotate("text", x=110, y=2.3, label = "182", size=2)+
  geom_segment(aes(x =  1.74665, y = 1.2, xend = 30, yend = 1.2), lty = "solid", lwd = 0.2, color = "black", #model
               arrow = grid::arrow(length = unit(0.1,"cm")))+
  geom_segment(aes(x = 141.994, y = 1.2, xend = 36, yend = 1.2), lty = "solid", lwd = 0.2, color = "black", #model
               arrow = grid::arrow(length = unit(0.1,"cm")))+
  annotate("text", x=12, y=1.3, label = "32", size=2)+
  annotate("text", x=70, y=1.3, label = "108", size=2)+
  geom_errorbarh(aes(xmax = area_weighted_flux + upper_95, xmin = area_weighted_flux - lower_95, height = .2, color = `Variability Source`))+
  geom_point(pch = 21, size = 6, color = "black")+
  labs(x = expression(paste("Area Corrected Global CH"[4]," Flux (g CH"[4]," m"^-2," yr"^-1,")")), y = "Emission Estimate Source")
  
p

ggsave("./Global_rate_comparison.tiff", device = "tiff", dpi = 1000, width = 20, height = 16, units = "cm")
  geom_errorbarh(aes(xmax = area_weighted_flux + upper_95, xmin = area_weighted_flux - lower_95, height = .2, color = group))+
  geom_point(pch = 21, size = 6, color = "black")+
  ylab("Emission Estaimte Source")+
  xlab("Area Corrected Global Methane Flux (g CH4 m-2 yr-1)")
  