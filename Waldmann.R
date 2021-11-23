
############PACKAGES

library(ggplot2)
library(ggpmisc)
library(gridExtra)
library(readxl)
library(Rmisc)
library(hrbrthemes)
library(GGally)
library(viridis)
library(ggplot2)
library(gridExtra)
library(readxl)
library(ggpubr)
library(Rmisc)
library(tidyverse)
library(plyr)

############DATA
library(readxl)
all <- read_excel("Pathology Reports/Waldmann/MHL 212448 Waldmann.xlsx", 
                   sheet = "PLOT")
data <- all
left <- all[ which(all$Side =='left'), ]
right <- all[ which(all$Side =='right'), ]
spleen <- all[ which(all$Side =='spleen'), ]
data <- all[ which(aall$Side =='right' & 'left'), ]
############PLOT

# Plot
Spleen.2 <- ggparcoord(spleen, columns = 8:20, groupColumn = 3, scale="std", showPoints = TRUE, alphaLines = .8, splineFactor = F, title = "Spleen") +
  #scale_color_viridis(discrete=TRUE) +
  theme_ipsum()+
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  #theme(axis.text.x=element_blank()) +
  #theme(axis.text.y=element_blank()) +
  theme(
    plot.title = element_text(size=13)) +
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank())
## scale = std globalminmax uniminmax center centerObs robust



tiff("Coord.Plots.tiff", units="in", width=28, height=14, res=250)
grid.arrange(Tumor.L.2, Tumor.L.4, Tumor.R.2, Tumor.R.4, Spleen.2, spleen.4, ncol = 2, nrow = 3)
dev.off()


tiff("4.groups.tiff", units="in", width=14, height=14, res=300)
grid.arrange(left.4, right.4, spleen.4, ncol = 1, nrow = 3)
dev.off()





########
my_mean = aggregate(data$'CD8/IFN-g per mm^2', by=list(data$Group), mean, na.rm=TRUE) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'CD8/IFN-g per mm^2' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

cd8 <- ggplot(data) + 
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Tumor: CD8/IFN-g per mm^2") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  #geom_jitter(aes(x = Group, y = PLT, color = Group), width = 0.1, show.legend=F)+
  theme_bw(base_size = 14) +
  geom_jitter(aes(x = Group, y = data$'CD8/IFN-g per mm^2', color = Side), width = 0.2, size = 3) +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x=element_text(angle=40,hjust=1)) 

cd8 <- cd8 + theme(legend.position="none")


tiff("Plot111.tiff", units="in", width=20, height=12, res=300)
grid.arrange(cd206, cd68, m1m2, cd8, cd11c, ncol = 3, nrow = 2)
dev.off()

tiff("cd8.tiff", units="in", width=13, height=10, res=150)
grid.arrange(cd8, ncol = 1, nrow = 1)
dev.off()


my_mean = aggregate(data$'3/4/21 BM Grade', by=list(data$'Group')) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'3/4/21 BM Grade', by=list(data$'Group'), FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean, my_CI, by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

ggplot(data) + 
  geom_jitter(aes(x = data$'Group', y = data$'3/4/21 BM Grade', color = data$'Group'), width = 0.2, size = 4) +
  geom_point(data = my_info, aes(x = my_info$'Group', y = my_info$mean), color = "grey", size = 5) +
  scale_y_continuous(name = "1st BM") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=2) +
  theme_bw(base_size = 18) +
  theme(axis.text.x=element_text(angle=15,hjust=1))+
  theme(axis.title.x=element_blank(), text = element_text(size = 20))#, legend.position="none")

