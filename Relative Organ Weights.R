

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
my_mean = aggregate(data$'Brain % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Brain % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Brain plot
Brain <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$"Brain % BW", color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Brain (as % of BW)")+#, limits=c(1.8, 2.3)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))


### Heart Weight
my_mean = aggregate(data$'Heart % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Heart % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Heart plot
Heart <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Heart % BW', color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Heart (as % of BW)")+#, limits=c(0.37, 0.56)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))

### Liver Weight
my_mean = aggregate(data$'Liver % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Liver % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Liver plot
Liver <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Liver % BW', color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Liver (as % of BW)")+#, limits=c(4, 5.3)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))

### Lung Weight
my_mean = aggregate(data$'Lung % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Lung % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Lung plot
Lung <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Lung % BW', color = Sex), width = 0.1, show.legend=F)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Lung (as % of BW)")+#, limits=c(0.5, 1.5)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))

### Spleen Weight
my_mean = aggregate(data$'Spleen % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Spleen % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Spleen plot
Spleen <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Spleen % BW', color = Sex), width = 0.1, show.legend=F) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Spleen (as % of BW)")+#, limits=c(0.04, 0.2)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))

### Kidney Weight
my_mean = aggregate(data$'Kidney % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Kidney % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Kidney plot
Kidney <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Kidney % BW', color = Sex), width = 0.1, show.legend=F) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Kidney (as % of BW)")+#, limits=c(0.8, 1.4)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))


### Graft Weight
my_mean = aggregate(data$'Allograft Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Allograft Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Graft plot
Graft <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Allograft Weight', color = Sex), width = 0.1, show.legend=F) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Graft Weight")+#, limits=c(0, 10)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=0.5) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0,hjust=0.5)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))


### Make multiple plots
tiff("BW.tiff", units="in", width=6, height=2.5, res=600)
grid.arrange(BW, 
             ncol = 2, nrow = 1)

dev.off()



### Make multiple plots
tiff("Relative Organ Weights.tiff", units="in", width=9, height=5, res=600)
grid.arrange(Brain, Spleen, Liver, Lung, Kidney, Heart, ncol = 3, nrow = 2)
dev.off()

#tiff("Relative Organ Weights.tiff", units="in", width=9, height=7.5, res=600)
#grid.arrange(BW, Graft, Brain, Spleen, Liver, Lung, Kidney, Heart, ncol = 3, nrow = 3)
#dev.off()

