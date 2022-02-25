
data <- read_excel("MHL 19-331-114 Efficacy.xlsx", sheet = "Chemistry")

library(ggplot2)
library(gridExtra)
library(readxl)
library(ggpubr)
library(readxl)

###Generate Data
data <- read_excel("ADME Tox 202.xlsx", sheet = "Chemistry")
CBC <- read_excel("ADME Tox 202.xlsx", sheet = "CBC")
AData <- read_excel("ADME Tox 202.xlsx", sheet = "Animal Data")

### ALB  
my_mean = aggregate(data$ALB, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$ALB , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(2.77)
my_info$ref.hi = c(4.8)

### ALB Plot
ALB <- ggplot(data) + 
  scale_y_continuous(name = "Albumin") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = ALB, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())


### ALP
my_mean = aggregate(data$ALP, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$ALP , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(41)
my_info$ref.hi = c(190)

### ALP plot
ALP <- ggplot(data) + 
  scale_y_continuous(name = "ALP") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = ALP, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())


### ALT
my_mean = aggregate(data$ALT, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$ALT , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(21)
my_info$ref.hi = c(66)

### ALT plot
ALT <- ggplot(data) + 
  scale_y_continuous(name = "ALT") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = ALT, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### BUN 
my_mean = aggregate(data$BUN, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$BUN , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(10.7)
my_info$ref.hi = c(34.2)

### BUN plot
BUN <- ggplot(data) + 
  scale_y_continuous(name = "BUN") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = BUN, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### CREATININE
my_mean = aggregate(data$CRE, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$CRE , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(0.2)
my_info$ref.hi = c(0.5)

### CREATININE plot
CRE <- ggplot(data) + 
  scale_y_continuous(name = "Creatinine") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = CRE, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### GLU 
my_mean = aggregate(data$GLU, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$GLU , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(115)
my_info$ref.hi = c(292)

### GLU plot
GLU <- ggplot(data) + 
  scale_y_continuous(name = "GLU") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = GLU, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### TP
my_mean = aggregate(data$TP, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$TP , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(4.6)
my_info$ref.hi = c(6.6)

### TP plot
TP <- ggplot(data) + 
  scale_y_continuous(name = "Total Protein") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = TP, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### GLOB
my_mean = aggregate(data$GLOB, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$GLOB , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(0.5)
my_info$ref.hi = c(3.1)

### GLOB plot
GLOB <- ggplot(data) + 
  scale_y_continuous(name = "GLOB") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = GLOB, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### Ca
my_mean = aggregate(data$Ca, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$Ca , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(9.77)
my_info$ref.hi = c(12.2)

### Ca plot
Ca <- ggplot(data) + 
  scale_y_continuous(name = "Ca2+") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = Ca, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### PHOS
my_mean = aggregate(data$PHOS, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$PHOS , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(9.1)
my_info$ref.hi = c(13.1)

### PHOS plot
PHOS <- ggplot(data) + 
  scale_y_continuous(name = "PHOS") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = PHOS, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### NaPlus
my_mean = aggregate(data$NaPlus, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$NaPlus , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(124)
my_info$ref.hi = c(174)

### NaPlus plot
NaPlus <- ggplot(data) + 
  scale_y_continuous(name = "Na+") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = NaPlus, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())

### KPlus
my_mean = aggregate(data$KPlus, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$KPlus , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(6)
my_info$ref.hi = c(11.8)

### KPlus plot
KPlus <- ggplot(data) + 
  scale_y_continuous(name = "K+") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#f5f5f5", width = 0, size=10) +
  geom_jitter(aes(x = Group, y = KPlus, color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "#a9a9a9", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "#a9a9a9", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank())


tiff("Chem.tiff", units="in", width=10, height=7, res=600)
#Add A:G ratio
grid.arrange(TP, ALB, GLOB, ALP, ALT, BUN, CRE, GLU, Ca, PHOS, KPlus, NaPlus, ncol = 4, nrow = 3)
dev.off()