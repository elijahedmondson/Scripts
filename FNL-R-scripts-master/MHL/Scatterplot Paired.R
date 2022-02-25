library(GGally)
library(ggplot2)
library(tidyverse)
library(gapminder)

theme_set(theme_bw(12))
plot<-humanized %>%
ggplot(aes(humanized$"Strain",humanized$"CD45%")) +
  geom_point(aes(color=humanized$"Exclude"), size=3) +
  geom_line(aes(group = humanized$"paired"), color = "grey") +
  scale_y_continuous(name = "% huCD45 Cells") + 
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 20))


ggsave("huCD45.percent.png")


tiff("Plot.tiff", units="in", width=12, height=6, res=150)
grid.arrange(plot, ncol = 1, nrow = 1)
dev.off()

