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


library(ggplot2)
library(gridExtra)
library(readxl)
library(ggpubr)
library(Rmisc)
library(tidyverse)
library(plyr)
library(readxl)

data <- read_excel("C:/Users/edmondsonef/Desktop/MHL Johnson Tollip.xlsx", 
                   sheet = "PJ List")

my_mean = aggregate(data$'% Vacuolated Phenotype', by=list(data$'Groups'), mean, na.rm=TRUE) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'% Vacuolated Phenotype', by=list(data$'Groups'), FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean, my_CI, by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

nucleus.vac <- ggplot(data) + 
  geom_point(data = my_info, aes(x = my_info$'Group', y = my_info$mean), color = "grey", size = 5) +
  scale_y_continuous(name = "% Vacuolated Phenotype") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  geom_jitter(aes(x = data$'Groups', y = data$'% Vacuolated Phenotype', color = data$'Sex'), width = 0.2, height = 0.00001, size = 2) +
  theme_bw(base_size = 18) +
  #theme(axis.text.x=element_text(angle=25,hjust=1))+
  theme(axis.title.x=element_blank(), text = element_text(size = 20))+
  theme(axis.title.x=element_blank(), legend.position="none")
nucleus.vac



############SD
############SD funs(mean, sem=sd(.)/sqrt(length(.)))
############SD
############SD
var <- 'Nucleus: Area µm^2 Average'

my_mean=aggregate(data$'Nucleus: Area µm^2 Average', by=list(data$Groups), mean, na.rm=TRUE) ; colnames(my_mean)=c("Group" , "mean")
my_sd=aggregate(data$'Nucleus: Area µm^2 Average', by=list(data$Groups), sd, na.rm=TRUE) ; colnames(my_sd)=c("Group" , "sd")
my_info=merge(my_mean, my_sd, by.x=1 , by.y=1)
my_info$se <- my_info$sd / sqrt(cdata$N)

nucleus.area <- ggplot(data) + 
  geom_point(data = my_info, aes(x = Group , y = my_info$mean), color = "grey", size = 5) +
  scale_y_continuous(name = 'Nucleus: Area µm^2 Average') + #, limits = c(15, 50)) +
  geom_errorbar(data = my_info, aes(x = Group, y = sd, ymin = mean - sd, ymax = mean + sd), color = "grey", width = 0.2 , size=2) +
  theme_bw(base_size = 18) +
  geom_jitter(aes(x = data$Groups, y = data$'Nucleus: Area µm^2 Average', color = data$Sex), width = 0.2, size = 4) +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 20), legend.position="none")
nucleus.area
############SEM
############SEM
############SEM 
############SEM
my_mean=aggregate(data$'Nucleus: Area µm^2 Average', by=list(data$Groups), mean, na.rm=TRUE) ; colnames(my_mean)=c("Group" , "mean")
my_sd=aggregate(data$'Nucleus: Area µm^2 Average', by=list(data$Groups), sd, na.rm=TRUE) ; colnames(my_sd)=c("Group" , "sd")
my_info=merge(my_mean, my_sd, by.x=1 , by.y=1)
my_info$se <- my_info$sd / sqrt(cdata$N)

ggplot(data) + 
  geom_point(data = my_info, aes(x = Group , y = my_info$mean), color = "grey", size = 5) +
  scale_y_continuous(name = "Nucleus: Area µm^2 Average") + #, limits = c(15, 50)) +
  geom_errorbar(data = my_info, aes(x = Group, y = se, ymin = mean - se, ymax = mean + se), color = "grey", width = 0.2 , size=2) +
  theme_bw(base_size = 18) +
  geom_jitter(aes(x = data$Groups, y = data$'Nucleus: Area µm^2 Average', color = data$Sex), width = 0.2, size = 6) +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 20), legend.position="none")




tiff("CD117.tiff", units="in", width=20, height=6, res=300)
grid.arrange(nucleus.vac, nucleus.area, nucleus.length, nucleus.max, ncol = 2, nrow = 2)
dev.off()