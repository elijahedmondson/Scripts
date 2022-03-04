##QuPath cell measurement summaries

library(dplyr)
library(tidyr)

data <- read.csv("C:/Users/edmondsonef/Desktop/measurements.csv")
tibble::as_tibble(data)
glimpse(data)
count(data, Image)
grouped.data <- data %>% group_by(Image) %>% summarise(funs(mean, median, sd))

data <- read.csv("C:/Users/edmondsonef/Desktop/measurements.csv")
tibble::as_tibble(data)
glimpse(data)


count(data, Image)

gr.data <- dplyr::group_by(data, Image)
median <- dplyr::summarise_each(gr.data, funs(median))
df <- data.frame(median)
write.csv(df, file = "C:/Users/edmondsonef/Desktop/median.csv")


#####CELL

total <- rbind(dataframeA, dataframeB)

dplyr::summarise_each(gr.data, funs(median))
dplyr::summarise_each(gr.data, funs(sd))


dplyr::group_by(data, SLIDE)