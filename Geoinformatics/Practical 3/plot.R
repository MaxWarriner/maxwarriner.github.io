library(tidyverse)
dat <- read_csv('reflectance.csv')
dat$Band <- factor(dat$Band, levels = c("ultrablue", "blue", "green", "red", "NIR", "SWIR", "MIR"))

ggplot(data = dat, aes(x = Band, y = Value, group = Location, color = Location)) + 
  geom_line() + 
  theme_dark() + 
  cowplot::theme_cowplot() + 
  scale_color_manual(values = c("Cork Harbour" = "blue", "Mudflats" = "brown", "Runway" = "grey3", "Terminal Building" = "grey", "Green Field 1" = "green1", "Green Field 2" = "green4","Non-Green Field" =  "yellow3"))
