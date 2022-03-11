##QuPath cell measurement summaries

library(dplyr)
library(tidyr)

data.vac <- read_excel("C:/Users/edmondsonef/Desktop/MHL Johnson Tollip.xlsx",
                   sheet = "vacuolated")

data3 <- full_join(data1, data2, by = "Image ID")

data <- read.csv("C:/Users/edmondsonef/Desktop/measurements.csv")
tibble::as_tibble(data)
glimpse(data)
num.tumor <- count(data, Image)

df <- count(data, Image, wre = data$`Vacuolated Tumor`)
df <- spread(data, col = `Vacuolated Tumor`, into = 6)

data %>% group_by(Image) %>% tally(grade = data$`Vacuolated Tumor`)

grouped.data <- data %>% group_by(Image) %>% summarise(data, funs(mean))

data <- read.csv("C:/Users/edmondsonef/Desktop/measurements.csv")
tibble::as_tibble(data)
glimpse(data)


count(data, Image)

gr.data <- dplyr::group_by(data, Image)
median <- dplyr::summarise_each(gr.data, funs(median))
df <- data.frame(median)
write.csv(df, file = "C:/Users/edmondsonef/Desktop/m.csv")


#####CELL

total <- rbind(dataframeA, dataframeB)

dplyr::summarise_each(gr.data, funs(median))
dplyr::summarise_each(gr.data, funs(sd))


dplyr::group_by(data, SLIDE)