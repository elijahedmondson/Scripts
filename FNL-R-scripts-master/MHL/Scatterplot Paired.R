library(GGally)
library(ggplot2)
library(tidyverse)
library(gapminder)

theme_set(theme_bw(12))
plot<-data %>%
ggplot(aes(data$"Groups",data$"WBC (11/10)")) +
  geom_point(aes(color=data$"Day"), size=3) +
  geom_line(aes(group = data$"paired"), color = "grey") +
  scale_y_continuous(name = "WBC count") + 
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 20))


ggsave("huCD45.percent.png")

setwd("C:/Users/edmondsonef/Desktop")
tiff("Plot.tiff", units="in", width=12, height=6, res=150)
plot
dev.off()

