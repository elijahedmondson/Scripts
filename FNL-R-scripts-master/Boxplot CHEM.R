###File
#data <- read_excel("C:/Users/edmondsonef/Desktop/ADME Tox 189.xlsx", sheet = "Chemistry")

library(ggplot2)
library(gridExtra)
library(readxl)
library(ggpubr)

###Generate Data

### ALB  
my_mean = aggregate(data$ALB, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$ALB , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(2.6)
my_info$ref.hi = c(4.6)

### ALB Plot
ALB <- ggplot(data) + 
  scale_y_continuous(name = "Albumin") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$ALB, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")


### ALP
my_mean = aggregate(data$ALP, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$ALP , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(41)
my_info$ref.hi = c(140)

### ALP plot
ALP <- ggplot(data) + 
  scale_y_continuous(name = "ALP") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$ALP, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")


### ALT
my_mean = aggregate(data$ALT, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$ALT , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(20)
my_info$ref.hi = c(65)

### ALT plot
ALT <- ggplot(data) + 
  scale_y_continuous(name = "ALT") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$ALT, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")

### BUN 
my_mean = aggregate(data$BUN, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$BUN , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(10.7)
my_info$ref.hi = c(29.9)

### BUN plot
BUN <- ggplot(data) + 
  scale_y_continuous(name = "BUN") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$BUN, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")

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
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$CRE, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")

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
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$GLU, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")

### TP
my_mean = aggregate(data$TP, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$TP , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(4.6)
my_info$ref.hi = c(6.9)

### TP plot
TP <- ggplot(data) + 
  scale_y_continuous(name = "Total Protein") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$TP, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")

### GLOB
my_mean = aggregate(data$GLOB, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$GLOB , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(.5)
my_info$ref.hi = c(3.2)

### GLOB plot
GLOB <- ggplot(data) + 
  scale_y_continuous(name = "GLOB") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$GLOB, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")

### Ca
my_mean = aggregate(data$Ca, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$Ca , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(9.634)
my_info$ref.hi = c(12.2)

### Ca plot
Ca <- ggplot(data) + 
  scale_y_continuous(name = "Ca2+") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$Ca, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")

### PHOS
my_mean = aggregate(data$PHOS, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$PHOS , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(8.9)
my_info$ref.hi = c(13.1)

### PHOS plot
PHOS <- ggplot(data) + 
  scale_y_continuous(name = "PHOS") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$PHOS, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")

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
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$NaPlus, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")

### KPlus
my_mean = aggregate(data$KPlus, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$KPlus , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(6)
my_info$ref.hi = c(11.7)

### KPlus plot
KPlus <- ggplot(data) + 
  scale_y_continuous(name = "K+") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$KPlus, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")

### TBIL
my_mean = aggregate(data$TBIL, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$TBIL , by=list(data$Group, color = data$Sex) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(0.1)
my_info$ref.hi = c(0.9)

### TBIL plot
TBIL <- ggplot(data) + 
  scale_y_continuous(name = "TBIL") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=4) +
  geom_jitter(aes(x = Group, y = data$TBIL, color = data$Sex), width = 0.1)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.position="none")

### Generate Multiplots
#ggarrange(TP, GLOB, ALB, ALP, ALT, BUN, CRE, TBIL,
#         labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
#          ncol = 2, nrow = 4)

#ggarrange(GLU, Ca, PHOS, NaPlus, KPlus,
#          labels = c("A", "B", "C", "D", "E"),
#          ncol = 2, nrow = 4)

#ggarrange(TP, ALB, BUN, CRE, Ca, PHOS, NaPlus, KPlus,
#          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
#          ncol = 2, nrow = 4)

ggarrange(TP, ALB, GLOB, ALP, ALT, BUN, CRE, GLU, Ca, PHOS, TBIL, KPlus, NaPlus,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N"),
          ncol = 3, nrow = 5)




### FOR ONLY 1 OB


library(ggplot2)
library(gridExtra)
library(readxl)
library(ggpubr)

###Generate Data

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
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$ALB), width = 0.1)+
  theme(axis.title.x=element_blank())


### ALP
my_mean = aggregate(data$ALP, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$ALP , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(61)
my_info$ref.hi = c(190)

### ALP plot
ALP <- ggplot(data) + 
  scale_y_continuous(name = "ALP") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$ALP), width = 0.1)+
  theme(axis.title.x=element_blank())


### ALT
my_mean = aggregate(data$ALT, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$ALT , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(28)
my_info$ref.hi = c(132)

### ALT plot
ALT <- ggplot(data) + 
  scale_y_continuous(name = "ALT") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$ALT), width = 0.1)+
  theme(axis.title.x=element_blank())

### BUN 
my_mean = aggregate(data$BUN, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$BUN , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(17)
my_info$ref.hi = c(27.1)

### BUN plot
BUN <- ggplot(data) + 
  scale_y_continuous(name = "BUN") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$BUN), width = 0.1)+
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
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$CRE), width = 0.1)+
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
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$GLU), width = 0.1)+
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
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$TP), width = 0.1)+
  theme(axis.title.x=element_blank())

### GLOB
my_mean = aggregate(data$GLOB, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$GLOB , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(0)
my_info$ref.hi = c(0.6)

### GLOB plot
GLOB <- ggplot(data) + 
  scale_y_continuous(name = "GLOB") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$GLOB), width = 0.1)+
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
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$Ca), width = 0.1)+
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
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$PHOS), width = 0.1)+
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
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$NaPlus), width = 0.1)+
  theme(axis.title.x=element_blank())

### KPlus
my_mean = aggregate(data$KPlus, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$KPlus , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(4.6)
my_info$ref.hi = c(8)

### KPlus plot
KPlus <- ggplot(data) + 
  scale_y_continuous(name = "K+") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$KPlus), width = 0.1)+
  theme(axis.title.x=element_blank())

### TBIL
my_mean = aggregate(data$TBIL, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$TBIL , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)
my_info$ref.low = c(0.1)
my_info$ref.hi = c(0.9)

### TBIL plot
TBIL <- ggplot(data) + 
  scale_y_continuous(name = "TBIL") +
  geom_errorbar(data = my_info, aes(x = Group, ymin = ref.low, ymax = ref.hi), color = "#CCFFFF", width = 0, size=10) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 4) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$TBIL), width = 0.1)+
  theme(axis.title.x=element_blank())

### Generate Multiplots
#ggarrange(TP, GLOB, ALB, ALP, ALT, BUN, CRE, TBIL,
#         labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
#          ncol = 2, nrow = 4)

#ggarrange(GLU, Ca, PHOS, NaPlus, KPlus,
#          labels = c("A", "B", "C", "D", "E"),
#          ncol = 2, nrow = 4)

#ggarrange(TP, ALB, BUN, CRE, Ca, PHOS, NaPlus, KPlus,
#          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
#          ncol = 2, nrow = 4)

ggarrange(TP, ALB, GLOB, ALP, ALT, BUN, CRE, GLU, Ca, PHOS, TBIL, KPlus, NaPlus, 
          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"),
          ncol = 3, nrow = 4)
