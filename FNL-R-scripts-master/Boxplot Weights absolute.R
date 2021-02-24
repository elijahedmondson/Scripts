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
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 5) +
  scale_y_continuous(name = "Body Weight (grams, 95% CI)") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=2) +
  theme_bw(base_size = 18) +
  geom_jitter(aes(x = Group, y = Weight), width = 0.1, size = 4)+
  theme(axis.title.x=element_blank())


### Brain Weight
my_mean = aggregate(data$'Brain Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Brain Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Brain plot
Brain <- ggplot(data) + 
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Brain Weight (grams)") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Brain Weight'), width = 0.1)+
  theme(axis.title.x=element_blank())


### Heart Weight
my_mean = aggregate(data$'Heart Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Heart Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Heart plot
Heart <- ggplot(data) + 
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Heart Weight (grams)") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Heart Weight'), width = 0.1)+
  theme(axis.title.x=element_blank())

### Liver Weight
my_mean = aggregate(data$'Liver Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Liver Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Liver plot
Liver <- ggplot(data) + 
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Liver Weight (grams)") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Liver Weight'), width = 0.1)+
  theme(axis.title.x=element_blank())

### Lung Weight
my_mean = aggregate(data$'Lung Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Lung Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Lung plot
Lung <- ggplot(data) + 
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Lung Weight (grams)") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Lung Weight'), width = 0.1)+
  theme(axis.title.x=element_blank())

### Spleen Weight
my_mean = aggregate(data$'Spleen Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Spleen Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Spleen plot
Spleen <- ggplot(data) + 
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Spleen Weight (grams)") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Spleen Weight'), width = 0.1) +
  theme(axis.title.x=element_blank())

### Kidney Weight
my_mean = aggregate(data$'Kidney Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Kidney Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Kidney plot
Kidney <- ggplot(data) + 
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Kidney Weight (grams)") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Kidney Weight'), width = 0.1) +
  theme(axis.title.x=element_blank())

### Xenograft Weight
my_mean = aggregate(data$'Xenograft Weight', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Xenograft Weight' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Xenograft plot
Xenograft <- ggplot(data) + 
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Xenograft") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Xenograft Weight'), width = 0.1) +
  theme(axis.title.x=element_blank())


### Make multiple plots
ggarrange(BW, Xenograft, Heart, Liver, Brain, Lung, Spleen,
          labels = c("A", "B", "C", "D", "E", "F", "G"),
          ncol = 2, nrow = 4)

### Make multiple plots
ggarrange(BW, Heart, Liver, Lung, Kidney,
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 3)




### Body Weight  
my_mean = aggregate(data$"Mammary Weight", by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$"Mammary Weight", by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### BW Plot
BW <- ggplot(data) + 
  geom_point(data = my_info, aes(x = Group , y = mean), color = "grey", size = 5) +
  scale_y_continuous(name = "Mammary Weight (grams, 95% CI)") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=2) +
  theme_bw(base_size = 18) +
  geom_jitter(aes(x = Group, y = data$"Mammary Weight"), width = 0.1, size = 4)+
  theme(axis.title.x=element_blank())




#############WHEN ONLY 1 OB


###Generate Data

### Body Weight  
my_mean = aggregate(data$Weight, by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$Weight , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### BW Plot
BW <- ggplot(data) + 
  geom_point(data = my_mean, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Body Weight (95% CI)") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = Weight), width = 0.1)+
  theme(axis.title.x=element_blank())


### Brain Weight
my_mean = aggregate(data$'Brain % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Brain % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Brain plot
Brain <- ggplot(data) + 
  geom_point(data = my_mean, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Brain (as % of BW)", limits=c(0, 4)) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Brain % BW'), width = 0.1)+
  theme(axis.title.x=element_blank())


### Heart Weight
my_mean = aggregate(data$'Heart % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Heart % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Heart plot
Heart <- ggplot(data) + 
  geom_point(data = my_mean, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Heart (as % of BW)", limits=c(0, 2)) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Heart % BW'), width = 0.1)+
  theme(axis.title.x=element_blank())

### Liver Weight
my_mean = aggregate(data$'Liver % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Liver % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Liver plot
Liver <- ggplot(data) + 
  geom_point(data = my_mean, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Liver (as % of BW)") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Liver % BW'), width = 0.1)+
  theme(axis.title.x=element_blank())

### Lung Weight
my_mean = aggregate(data$'Lung % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Lung % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Lung plot
Lung <- ggplot(data) + 
  geom_point(data = my_mean, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Lung (as % of BW)") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Lung % BW'), width = 0.1)+
  theme(axis.title.x=element_blank())

### Spleen Weight
my_mean = aggregate(data$'Spleen % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Spleen % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Spleen plot
Spleen <- ggplot(data) + 
  geom_point(data = my_mean, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Spleen (as % of BW)", limits=c(0, 2)) +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Spleen % BW'), width = 0.1) +
  theme(axis.title.x=element_blank())

### Kidney Weight
my_mean = aggregate(data$'Kidney % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Kidney % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Kidney plot
Kidney <- ggplot(data) + 
  geom_point(data = my_mean, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Kidney (as % of BW)") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Kidney % BW'), width = 0.1) +
  theme(axis.title.x=element_blank())

### Xenograft Weight
my_mean = aggregate(data$'Xenograft % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Xenograft % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Xenograft plot
Xenograft <- ggplot(data) + 
  geom_point(data = my_mean, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Xenograft (as % of BW)") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Xenograft % BW'), width = 0.1) +
  theme(axis.title.x=element_blank())


### Make multiple plots
ggarrange(BW, Xenograft, Heart, Liver, Brain, Lung, Spleen,
          labels = c("A", "B", "C", "D", "E", "F", "G"),
          ncol = 2, nrow = 4)

### Make multiple plots
ggarrange(BW, Heart, Liver, Lung, Kidney,
          labels = c("A", "B", "C", "D", "E"),
          ncol = 2, nrow = 3)


### Pancreas Weight
my_mean = aggregate(data$'Pancreas % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Pancreas % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### Pancreas plot
Pancreas <- ggplot(data) + 
  geom_point(data = my_mean, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "Pancreas (as % of BW)") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'Pancreas % BW'), width = 0.1) +
  theme(axis.title.x=element_blank())


### LN, Mesenteric Weight
my_mean = aggregate(data$'LN, Mesenteric % BW', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'LN, Mesenteric % BW' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### LN, Mesenteric plot
LN.Mesenteric <- ggplot(data) + 
  geom_point(data = my_mean, aes(x = Group , y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "LN, Mesenteric (as % of BW)") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'LN, Mesenteric % BW'), width = 0.1) +
  theme(axis.title.x=element_blank())

### Make multiple plots
ggarrange(BW, Pancreas, Heart, Liver, Brain, Lung, Spleen, LN.Mesenteric,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 2, nrow = 4)

### OTHER#################################################################################################################################### 
my_mean = aggregate(data$'CD4: Percent Positive Nuclei', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'CD4: Percent Positive Nuclei' , by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("'Ext ID '" , "CI")
my_info = merge(my_mean , my_CI , by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

### BW Plot
ggplot(data) + 
  geom_point(data = my_info, aes(x = Group, y = mean), color = "grey", size = 3) +
  scale_y_continuous(name = "CD4: Percent Positive Nuclei") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'CD4: Percent Positive Nuclei'), width = 0.1)+
  theme(axis.title.x=element_blank())




#Calculation of mean and sd of each group ?
my_mean = aggregate(data$'CD4:CD8 Ratio', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_sd=aggregate(data$'CD4:CD8 Ratio' , by=list(data$Group) , sd) ; colnames(my_sd)=c("Group" , "sd")
my_info=merge(my_mean , my_sd , by.x=1 , by.y=1)

ggplot(data) + 
  geom_point(data = my_info, aes(x = Group, y = mean), color = "grey", size = 6) +
  scale_y_continuous(name = "CD4:CD8 Ratio") +
  geom_errorbar(data = my_info, aes(x = Group,  y = sd, ymin = mean - sd, ymax = mean + sd), color = "grey", width = 0.4 , size=1) +
  theme_bw() +
  geom_jitter(aes(x = Group, y = data$'CD4:CD8 Ratio'), width = 0.2, size = 4)+
  theme(axis.title.x=element_blank())


