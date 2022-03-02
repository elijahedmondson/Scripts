##QuPath cell measurement summaries

library(dplyr)
library(tidyr)

data <- read.csv("C:/Users/edmondsonef/Desktop/measurements.csv")
tibble::as_tibble(data)
glimpse(data)
count(data, Image)
grouped.data <- data %>% group_by(Image) %>% summarise(funs(mean, median, sd))

data <- read.csv("C:/Users/edmondsonef/Desktop/cells.csv")
tibble::as_tibble(data)
glimpse(data)


count(data, Image)

gr.data <- dplyr::group_by(data, Image)
sd <- dplyr::summarise_each(gr.data, funs(sd))
df <- data.frame(sd)
write.csv(df, file = "C:/Users/edmondsonef/Desktop/sd.csv")

dplyr::summarise_each(gr.data, funs(median))
dplyr::summarise_each(gr.data, funs(sd))


dplyr::group_by(data, SLIDE)