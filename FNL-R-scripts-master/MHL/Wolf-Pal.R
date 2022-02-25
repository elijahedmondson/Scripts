
library(readxl)
data <- read_excel("C:/Users/edmondsonef/Desktop/MHL 220584A Pal Wolf.xlsx", sheet = "Summary")


library(ggplot2)
library(gridExtra)
library(readxl)
library(ggpubr)
library(Rmisc)
library(tidyverse)
library(plyr)

my_mean = aggregate(data$'Margin CD86 Cells per ??m²', by=list(data$'Group'), mean, na.rm=TRUE) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Margin CD86 Cells per ??m²', by=list(data$'Group'), FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean, my_CI, by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

Image1 <- ggplot(data) + 
  geom_point(data = my_info, aes(x = my_info$'Group', y = my_info$mean), color = "grey", size = 5) +
  scale_y_continuous(name = "Margin CD86: Cells per ??m²") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  geom_jitter(aes(x = data$'Group', y = data$'Margin CD86 Cells per ??m²', color = data$'Day', shape = data$'Batch'), width = 0.2, height = 0.1, size = 2) +
  theme_bw(base_size = 18) +
  theme(axis.text.x=element_text(angle=35,hjust=1))+
  theme(axis.title.x=element_blank(), text = element_text(size = 20))#, legend.position="none")

tiff("CD86.tiff", units="in", width=9, height=6, res=300)
grid.arrange(Image1, ncol = 1, nrow = 1)
dev.off()


my_mean = aggregate(data$'Margin F4/80 Cells per ??m²', by=list(data$'Group'), mean, na.rm=TRUE) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Margin F4/80 Cells per ??m²', by=list(data$'Group'), FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean, my_CI, by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

Image2 <- ggplot(data) + 
  geom_point(data = my_info, aes(x = my_info$'Group', y = my_info$mean), color = "grey", size = 5) +
  scale_y_continuous(name = "Margin F4/80: Cells per ??m²") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  geom_jitter(aes(x = data$'Group', y = data$'Margin F4/80 Cells per ??m²', color = data$'Day', shape = data$'Batch'), width = 0.2, height = 0.1, size = 2) +
  theme_bw(base_size = 18) +
  theme(axis.text.x=element_text(angle=35,hjust=1))+
  theme(axis.title.x=element_blank(), text = element_text(size = 20))#, legend.position="none")
  
tiff("F480.tiff", units="in", width=9, height=6, res=300)
grid.arrange(Image2, ncol = 1, nrow = 1)
dev.off()



my_mean = aggregate(data$'Margin Ly6G Cells per ??m²', by=list(data$'Group'), mean, na.rm=TRUE) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Margin Ly6G Cells per ??m²', by=list(data$'Group'), FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean, my_CI, by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

Image3 <- ggplot(data) + 
  geom_point(data = my_info, aes(x = my_info$'Group', y = my_info$mean), color = "grey", size = 5) +
  scale_y_continuous(name = "Margin Ly6G: Cells per ??m²") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  geom_jitter(aes(x = data$'Group', y = data$'Margin Ly6G Cells per ??m²', color = data$'Day', shape = data$'Batch'), width = 0.2, height = 0.1, size = 2) +
  theme_bw(base_size = 18) +
  theme(axis.text.x=element_text(angle=35,hjust=1))+
  theme(axis.title.x=element_blank(), text = element_text(size = 20))#, legend.position="none")
  
  tiff("Ly6G.tiff", units="in", width=9, height=6, res=300)
grid.arrange(Image3, ncol = 1, nrow = 1)
dev.off()



my_mean = aggregate(data$'Margin CD86/F480 Cells per ??m²', by=list(data$'Group'), mean, na.rm=TRUE) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Margin CD86/F480 Cells per ??m²', by=list(data$'Group'), FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean, my_CI, by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

Image4 <- ggplot(data) + 
  geom_point(data = my_info, aes(x = my_info$'Group', y = my_info$mean), color = "grey", size = 5) +
  scale_y_continuous(name = "Margin CD86/F480: Cells per ??m²") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  geom_jitter(aes(x = data$'Group', y = data$'Margin CD86/F480 Cells per ??m²', color = data$'Day', shape = data$'Batch'), width = 0.2, height = 0.1, size = 2) +
  theme_bw(base_size = 18) +
  theme(axis.text.x=element_text(angle=35,hjust=1))+
  theme(axis.title.x=element_blank(), text = element_text(size = 20))#, legend.position="none")
  
tiff("CD86+F480.tiff", units="in", width=9, height=6, res=300)
grid.arrange(Image4, ncol = 1, nrow = 1)
dev.off()






my_mean = aggregate(data$'Margin Total Cells per ??m²', by=list(data$'Group'), mean, na.rm=TRUE) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Margin Total Cells per ??m²', by=list(data$'Group'), FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean, my_CI, by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

Image5 <- ggplot(data) + 
  geom_point(data = my_info, aes(x = my_info$'Group', y = my_info$mean), color = "grey", size = 5) +
  scale_y_continuous(name = "Margin Total: Cells per ??m²") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  geom_jitter(aes(x = data$'Group', y = data$'Margin Total Cells per ??m²', color = data$'Day', shape = data$'Batch'), width = 0.2, height = 0.1, size = 2) +
  theme_bw(base_size = 18) +
  theme(axis.text.x=element_text(angle=35,hjust=1))+
  theme(axis.title.x=element_blank(), text = element_text(size = 20))#, legend.position="none")

tiff("Total.tiff", units="in", width=9, height=6, res=300)
grid.arrange(Image5, ncol = 1, nrow = 1)
dev.off()


