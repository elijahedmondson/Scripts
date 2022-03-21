
#data <- read_excel("MHL 19-331-114 Efficacy.xlsx", sheet = "CBC")

library(ggplot2)
library(gridExtra)
library(readxl)
#library(ggpubr)

###Generate Data

### WBC  
my_mean = aggregate(data$WBC, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$WBC , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(2.71)
my_info$ref.hi = c(12.33)

### WBC Plot
WBC <- ggplot(data) + 
  scale_y_continuous(name = "White Blood Cells") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = WBC, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())


### Neutrophils
my_mean = aggregate(data$NE, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$NE , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(0.46)
my_info$ref.hi = c(6.04)

### Nuetrophils plot
NE <- ggplot(data) + 
  scale_y_continuous(name = "Neutrophils") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = NE, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())


### Lymphocyes
my_mean = aggregate(data$LY, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$LY , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(1.25)
my_info$ref.hi = c(6.99)

### Lymphocyes plot
LY <- ggplot(data) + 
  scale_y_continuous(name = "Lymphocyes") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = LY, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### Monocytes 
my_mean = aggregate(data$MO, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$MO , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(0)
my_info$ref.hi = c(1.04)

### Monocytes plot
MO <- ggplot(data) + 
  scale_y_continuous(name = "Monocytes") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = MO, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### Eosinophils
my_mean = aggregate(data$EO, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$EO , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(0)
my_info$ref.hi = c(0.25)

### Eosinophils plot
EO <- ggplot(data) + 
  scale_y_continuous(name = "Eosinophils") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = EO, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### Basophils 
my_mean = aggregate(data$BA, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$BA , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(0)
my_info$ref.hi = c(0.23)

### Basophils plot
BA <- ggplot(data) + 
  scale_y_continuous(name = "Basophils") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = BA, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### Hematocrit
my_mean = aggregate(data$HCT, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$HCT , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(36.85)
my_info$ref.hi = c(50.71)

### Hematocrit plot
HCT <- ggplot(data) + 
  scale_y_continuous(name = "Hematocrit")+#, limits=c(30, 60)) +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = HCT, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### Platelets
my_mean = aggregate(data$PLT, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$PLT , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(305.79)
my_info$ref.hi = c(1444.39)

### Platelets plot
PLT <- ggplot(data) + 
  scale_y_continuous(name = "Platelets") + 
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = PLT, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())



### Generate Multiplots

tiff("WBC.tiff", units="in", width=8, height=4.667, res=600)
grid.arrange(WBC, NE, LY, MO, EO, BA, ncol = 3, nrow = 2)
dev.off()


## RBC ###
## RBC ###
## RBC ###
## RBC ###
## RBC ###
## RBC ###
## RBC ###
## RBC ###
## RBC ###
## RBC ###
## RBC ###
## RBC ###
## RBC ###

### Hematocrit
my_mean = aggregate(data$HCT, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$HCT , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(36.85)
my_info$ref.hi = c(50.71)

### Hematocrit plot
HCT <- ggplot(data) + 
  scale_y_continuous(name = "Hematocrit")+#, limits=c(30, 60)) +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = HCT, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())


### Reticulocytes
my_mean = aggregate(data$Retics, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$Retics , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(6)
my_info$ref.hi = c(6)

### Reticulocytes plot
Retics <- ggplot(data) + 
  scale_y_continuous(name = "Reticulocytes") + 
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = Retics, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### MCH
my_mean = aggregate(data$MCH, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$MCH , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(14)
my_info$ref.hi = c(20)

### MCH plot
MCH <- ggplot(data) + 
  scale_y_continuous(name = "MCH")+#, limits=c(12, 21)) +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = MCH, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### MCV
my_mean = aggregate(data$MCV, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$MCV , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(47)
my_info$ref.hi = c(67)

### MCV plot
MCV <- ggplot(data) + 
  scale_y_continuous(name = "MCV")+#, limits=c(40, 68)) + 
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = MCV, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### MCHC
my_mean = aggregate(data$MCHC, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$MCHC, by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(29)
my_info$ref.hi = c(33)

### MCHC plot
MCHC <- ggplot(data) + 
  scale_y_continuous(name = "MCHC")+#, limits=c(28, 34)) +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = MCHC, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### Hb
my_mean = aggregate(data$Hb, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$Hb, by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(11)
my_info$ref.hi = c(15.6)

### Hb plot
Hb <- ggplot(data) + 
  scale_y_continuous(name = "Hb")+#, limits=c(10, 16)) +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = Hb, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### RBC
my_mean = aggregate(data$RBC, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$RBC, by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(6)
my_info$ref.hi = c(11)

### RBC plot
RBC <- ggplot(data) + 
  scale_y_continuous(name = "RBC")+#, limits=c(5.5, 12)) +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = RBC, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### RDW
my_mean = aggregate(data$RDW, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$RDW, by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(13)
my_info$ref.hi = c(20)

### RDW plot
RDW <- ggplot(data) + 
  scale_y_continuous(name = "RDW")+#, limits=c(12, 21)) +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = RDW, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

tiff("RBC.tiff", units="in", width=8, height=7, res=600)
grid.arrange(HCT, Hb, MCH, MCHC, MCV, RBC, RDW, Retics, ncol = 3, nrow = 3)
dev.off()

###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS
###PLATELETS

### Platelets
my_mean = aggregate(data$PLT, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$PLT , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(305.79)
my_info$ref.hi = c(1444.39)

### Platelets plot
PLT <- ggplot(data) + 
  scale_y_continuous(name = "Platelets") + 
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = PLT, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())


### MPV
my_mean = aggregate(data$MPV, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$MPV , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(4)
my_info$ref.hi = c(10.6)

### Platelets plot
MPV <- ggplot(data) + 
  scale_y_continuous(name = "MPV") + 
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = MPV, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.3 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

tiff("Platelets.tiff", units="in", width=8, height=2.33, res=600)
grid.arrange(PLT, MPV, ncol = 3, nrow = 1)
dev.off()
