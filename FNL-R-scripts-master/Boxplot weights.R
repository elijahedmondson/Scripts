

library(ggplot2)
library(gridExtra)
library(readxl)
library(ggpubr)



###Generate Data

### Body Weight  
my_mean = aggregate(data$Weight, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$Weight , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### BW Plot
BW <- ggplot(data) + 
  scale_y_continuous(name = "Body Weight (95% CI)") +
  theme_bw() +
  geom_jitter(aes(x = Group, y = Weight, color = data$Sex), width = 0.1, size = 2)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2.5) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  #theme(axis.text.x=element_text(angle=0,hjust=0)) +
  #theme(axis.title.x=element_blank(), text = element_text(size = 8), legend.title = element_blank(), legend.position = c(0.8, 0.2), legend.background = element_rect(linetype= "solid", size = 0.1, colour = "black"))
  theme(axis.title.x=element_blank(), text = element_text(size = 12), legend.position="none")

### Brain Weight
my_mean = aggregate(data$'Brain % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Brain % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Brain plot
Brain <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Brain % BW', color = data$Sex), width = 0.1, size = 2)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2.5) +
  scale_y_continuous(name = "Brain (as % of BW)") + #, limits=c(0, 4)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  #theme(axis.text.x=element_text(angle=0,hjust=0)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 12), legend.position="none")


### Heart Weight
my_mean = aggregate(data$'Heart % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Heart % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Heart plot
Heart <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Heart % BW', color = data$Sex), width = 0.1, size = 2)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2.5) +
  scale_y_continuous(name = "Heart (as % of BW)") + #, limits=c(0, 2)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  #theme(axis.text.x=element_text(angle=0,hjust=0)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 12), legend.position="none")

### Liver Weight
my_mean = aggregate(data$'Liver % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Liver % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Liver plot
Liver <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Liver % BW', color = data$Sex), width = 0.1, size = 2)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2.5) +
  scale_y_continuous(name = "Liver (as % of BW)") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  #theme(axis.text.x=element_text(angle=0,hjust=0)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 12), legend.position="none")

### Lung Weight
my_mean = aggregate(data$'Lung % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Lung % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Lung plot
Lung <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Lung % BW', color = data$Sex), width = 0.1, size = 2)+
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2.5) +
  scale_y_continuous(name = "Lung (as % of BW)") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  #theme(axis.text.x=element_text(angle=0,hjust=0)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 12), legend.position="none")

### Spleen Weight
my_mean = aggregate(data$'Spleen % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Spleen % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Spleen plot
Spleen <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Spleen % BW', color = data$Sex), width = 0.1, size = 2) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2.5) +
  scale_y_continuous(name = "Spleen (as % of BW)") + #, limits=c(0, 2)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  #theme(axis.text.x=element_text(angle=0,hjust=0)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 12), legend.position="none")

### Kidney Weight
my_mean = aggregate(data$'Kidney % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Kidney % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Kidney plot
Kidney <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Kidney % BW', color = data$Sex), width = 0.1, size = 2) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2.5) +
  scale_y_continuous(name = "Kidney (as % of BW)") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  #theme(axis.text.x=element_text(angle=0,hjust=0)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 12), legend.position="none")

### Xenograft Weight
my_mean = aggregate(data$'Xenograft Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Xenograft Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Xenograft Weight plot
Xenograft <- ggplot(data) + 
  geom_jitter(aes(x = Group, y = data$'Xenograft Weight', color = data$Sex), width = 0.1, size = 2) +
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2.5) +
  scale_y_continuous(name = "Xenograft Weight") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  #theme(axis.text.x=element_text(angle=0,hjust=0)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 12), legend.position="none")

### Make multiple plots
### Make multiple plots
### Make multiple plots
ggarrange(BW, Xenograft, Brain, Spleen, Liver, Lung, Kidney, Heart,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 2, nrow = 4)


### Make multiple plots
### Make multiple plots
### Make multiple plots
ggarrange(BW, Brain, Spleen, Liver, Lung, Kidney, Heart,
          labels = c("A", "B", "C", "D", "E", "F", "G"),
          ncol = 2, nrow = 4)

### Make multiple plots
### Make multiple plots
### Make multiple plots
ggarrange(BW, Spleen, Liver, Lung, Kidney, Heart,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)

### Allograft Weight
my_mean = aggregate(data$'Allograft Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Allograft Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Allograft plot
Allograft <- ggplot(data) + 
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 2) +
  scale_y_continuous(name = "Allograft Weight") + #, limits=c(0, 10)) +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Allograft Weight'), width = 0.1, size = 1) +
  theme(axis.text.x=element_text(angle=0,hjust=0)) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8))

### Make multiple plots
### Make multiple plots
### Make multiple plots
ggarrange(BW, Allograft, Brain, Spleen, Liver, Lung, Kidney, Heart,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 2, nrow = 4)




### Make multiple plots
ggarrange(BW, Spleen, Liver, Lung, Kidney, Heart,
          labels = c("A", "B", "C", "D", "E", "F"),
          ncol = 2, nrow = 3)












######


library(ggplot2)
library(gridExtra)
library(readxl)
library(ggpubr)



###Generate Data

### Body Weight  
my_mean = aggregate(data$Weight, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")

### BW Plot
BW <- ggplot(data) + 
  scale_y_continuous(name = "Body Weight (95% CI)") +
  theme_bw() +
  geom_jitter(aes(x = Group, y = Weight), width = 0.1, size = 4)+
  theme(axis.title.x=element_blank(), text = element_text(size = 8))


### Brain Weight
my_mean = aggregate(data$'Brain % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")

### Brain plot
Brain <- ggplot(data) + 
  scale_y_continuous(name = "Brain (as % of BW)") + #, limits=c(0, 4)) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Brain % BW'), width = 0.1, size = 1)+
  theme(axis.title.x=element_blank(), text = element_text(size = 8))


### Heart Weight
my_mean = aggregate(data$'Heart % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")

### Heart plot
Heart <- ggplot(data) + 
  scale_y_continuous(name = "Heart (as % of BW)") + #, limits=c(0, 2)) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Heart % BW'), width = 0.1, size = 1)+
  theme(axis.title.x=element_blank(), text = element_text(size = 8))

### Liver Weight
my_mean = aggregate(data$'Liver % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")

### Liver plot
Liver <- ggplot(data) + 
  scale_y_continuous(name = "Liver (as % of BW)") +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Liver % BW'), width = 0.1, size = 1)+
  theme(axis.title.x=element_blank(), text = element_text(size = 8))

### Lung Weight
my_mean = aggregate(data$'Lung % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")

### Lung plot
Lung <- ggplot(data) + 
  scale_y_continuous(name = "Lung (as % of BW)") +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Lung % BW'), width = 0.1, size = 1)+
  theme(axis.title.x=element_blank(), text = element_text(size = 8))

### Spleen Weight
my_mean = aggregate(data$'Spleen % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")

### Spleen plot
Spleen <- ggplot(data) + 
  scale_y_continuous(name = "Spleen (as % of BW)") + #, limits=c(0, 2)) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Spleen % BW'), width = 0.1, size = 1) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8))

### Kidney Weight
my_mean = aggregate(data$'Kidney % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")

### Kidney plot
Kidney <- ggplot(data) + 
  scale_y_continuous(name = "Kidney (as % of BW)") +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Kidney % BW'), width = 0.1, size = 1) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8))

### Mesenteric LN Weight
my_mean = aggregate(data$'Mesenteric LN % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")

### Mesenteric LN  plot
Mesenteric <- ggplot(data) + 
  scale_y_continuous(name = "Mesenteric LN % BW") + #, limits=c(0, 10)) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Mesenteric LN % BW'), width = 0.1, size = 1) +
  theme(axis.title.x=element_blank(), text = element_text(size = 8))


### Make multiple plots
ggarrange(BW, Mesenteric, Brain, Spleen, Liver, Lung, Kidney, Heart,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 2, nrow = 4)




