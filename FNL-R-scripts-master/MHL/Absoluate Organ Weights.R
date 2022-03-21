
#data <- read_excel("C:/Users/edmondsonef/Desktop/MHL 19-331-121 Efficacy.xlsx", sheet = "19-331-121 Chem")
library(ggplot2)
library(gridExtra)
library(readxl)



### Body Weight  
my_mean = aggregate(data$Weight, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$Weight , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### BW Plot
BW <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = Weight, color = Sex), width = 0.1, size = 2, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Body Weight (95% CI)")+#, limits=c(21, 26.5)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))


### Brain Weight
my_mean = aggregate(data$'Brain Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Brain Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Brain plot
Brain <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$"Brain Weight", color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Brain Weight (grams)")+#, limits=c(0.4, 0.55)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))


### Heart Weight
my_mean = aggregate(data$'Heart Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Heart Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Heart plot
Heart <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Heart Weight', color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Heart Weight (grams)")+#, limits=c(0.06, 0.15)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))

### Liver Weight
my_mean = aggregate(data$'Liver Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Liver Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Liver plot
Liver <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Liver Weight', color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Liver Weight (grams)")+#, limits=c(0.9, 1.3)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))

### Lung Weight
my_mean = aggregate(data$'Lung Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Lung Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Lung plot
Lung <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Lung Weight', color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Lung Weight (grams)")+#, limits=c(0.1, 0.36)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))

### Spleen Weight
my_mean = aggregate(data$'Spleen Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Spleen Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Spleen plot
Spleen <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Spleen Weight', color = Sex), width = 0.1, show.legend=F) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Spleen Weight (grams)")+#, limits=c(0.009, 0.1)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))

### Kidney Weight
my_mean = aggregate(data$'Kidney Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Kidney Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Kidney plot
Kidney <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Kidney Weight', color = Sex), width = 0.1, show.legend=F) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Kidney Weight (grams)")+#, limits=c(0.18, 0.33)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))


### Graft Weight
my_mean = aggregate(data$'Allograft Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Allograft' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Graft plot
Graft <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Allograft Weight', color = Sex), width = 0.1, show.legend=F) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Graft")+#, limits=c(0, 10)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))



### Make multiple plots
tiff("Absolute Organ Weights.tiff", units="in", width=9, height=5, res=600)
grid.arrange(Brain, Spleen, Liver, Lung, Kidney, Heart,
             ncol = 3, nrow = 2)
dev.off()

#tiff("Absolute Organ Weights.tiff", units="in", width=9, height=7.5, res=600)
#grid.arrange(BW, Graft, Brain, Spleen, Liver, Lung, Kidney, Heart, ncol = 3, nrow = 3)
#dev.off()