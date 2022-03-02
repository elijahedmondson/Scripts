##QuPath cell measurement summaries


data <- read.csv("C:/Users/edmondsonef/Desktop/cells.csv")

Li
data$Nucleus..Area.µm.2

dplyr::summarise_each(data, funs(mean))
dplyr::summarise_each(data, funs(median))
dplyr::summarise_each(data, funs(sd))


dplyr::group_by(data, SLIDE)