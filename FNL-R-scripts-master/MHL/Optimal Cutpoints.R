library(OptimalCutpoints)
library(ggplot2)
library(gridExtra)
library(readxl)
library(ggpubr)
library(Rmisc)

summary.optimal.cutpoints(data$`CD8 (cells per mm^2)`)

optimal.cutpoints(data$`CD8 (cells per mm^2)`, methods = "ROC01")

optimal.cutpoints(data$`CD8 (cells per mm^2)`, status, tag.healthy, methods = "ROC01", data = data, direction = c("<", ">"), 
                  categorical.cov = NULL, pop.prev = NULL, control = control.cutpoints(), 
                  ci.fit = FALSE, conf.level = 0.95, trace = FALSE)


##LOWESS Smoothing

plot(data$`PD1 (cells per mm^2)`, data$dist, main="Locally Weighted Scaterpolot Smoothing: PD1+ cells per mm^2")
lines(lowess(data$`PD1 (cells per mm^2)`, data$dist), col=2)
lines(lowess(data$`PD1 (cells per mm^2)`, data$dist, f=.2), col=3)
legend(5, 620, c(paste("f=", c("2/3", ".2"))), lty=1, col=2:3)

####################
####################
####################
####################

my_mean = aggregate(data$'Inflammatory Infiltrates (average)', by=list(data$Group), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'Inflammatory Infiltrates (average)', by=list(data$Group) , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean, my_CI, by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

ggplot(data) + 
  geom_point(data = my_info, aes(x = Group, y = my_info$mean), color = "grey", size = 7) +
  scale_y_continuous(name = "Inflammatory Infiltrates (average)") +
  geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  theme_bw(base_size = 18) +
  geom_jitter(aes(x = data$Group, y = data$'Inflammatory Infiltrates (average)', color= data$'Kidney* PCR'), width = 0.2, size = 4) +
  #theme(axis.text.x=element_text(angle=25,hjust=1)) +
  theme(axis.title.x=element_blank(), legend.title = element_blank())

+
  geom_text( label=data$`Label`, size = 2,
             nudge_x = 0, nudge_y = 20, 
             check_overlap = T)



####################
####################
####################
####################

my_mean = aggregate(data$'PD1 (cells per mm^2)', by=list(data$'Tumor'), mean) ; colnames(my_mean)=c("Group" , "mean")
my_CI = aggregate(data$'PD1 (cells per mm^2)', by=list(data$'Tumor') , FUN = function(x) t.test(x)$conf.int) ; colnames(my_CI)=c("Group" , "CI")
my_info = merge(my_mean, my_CI, by.x=1 , by.y=1)
my_info$CIdiff = ((my_CI$CI[,2] - my_CI$CI[,1])/2)

ggplot(data, aes(x=data$'Tumor', y=data$'PD1 (cells per mm^2)', color=Group)) +
  #geom_point(shape=0, size = 4) +
  scale_y_continuous(name = "PD1 (cells per mm^2)") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  scale_colour_hue(l=50)+
  theme_bw(base_size = 22) +
  geom_jitter(aes(x = data$'Tumor', y = data$'PD1 (cells per mm^2)', color=Group), shape = 1, width = 0.1, size = 4) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(axis.title.x=element_blank(), legend.title = element_blank(), legend.justification=c("left", "top")) +
  #geom_label(label=data$`Label`, label.size = .01, nudge_x = 1, nudge_y = 5)
  geom_text(label=data$`Label`, hjust = 0, size = 1.7,
            nudge_x = .2, nudge_y = 0, 
             check_overlap = T)



png("PD1.png", width = 2000, height = 500, res = 300)
ggplot(data, aes(x=data$'Tumor', y=data$'PD1 (cells per mm^2)', color=Group)) +
  #geom_point(shape=0, size = 4) +
  scale_y_continuous(name = "PD1 (cells per mm^2)") +
  #geom_errorbar(data = my_info, aes(x = Group, y = CIdiff, ymin = mean - CIdiff, ymax = mean + CIdiff), color = "grey", width = 0.2 , size=1) +
  scale_colour_hue(l=50)+
  theme_bw(base_size = 24) +
  geom_jitter(aes(x = data$'Tumor', y = data$'PD1 (cells per mm^2)', color=Group), shape = 1, width = 0.1, size = 4) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  theme(axis.title.x=element_blank(), legend.title = element_blank(), legend.justification=c("left", "top")) +
  #geom_label(label=data$`Label`, label.size = .01, nudge_x = 1, nudge_y = 5)
  geom_text(label=data$`Label`, hjust = 0, size = 2,
            nudge_x = .1, nudge_y = 0, 
            check_overlap = T)
dev.off()


library(ggplot2)
library(ggpmisc)
library(ggplot2)
library(gridExtra)
library(readxl)
library(ggpubr)
library(Rmisc)


my.formula <- y ~ x
ggplot(data = data, aes(x = data$'CD8 (cells per mm^2)', y = data$'CD3 (cells per mm^2)', color = data$"Tumor Code"), na.rm=TRUE) +
  #geom_smooth(method = "lm", se=FALSE, color="red", formula = my.formula, na.rm=TRUE) +
  #stat_poly_eq(formula = my.formula, 
  #             aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
  #             parse = TRUE, na.rm=TRUE) +  
  geom_point(na.rm=TRUE, size = 3)+
  scale_y_continuous(name = "CD3 (cells per mm^2)") +
  scale_x_continuous(name = "CD8 (cells per mm^2)") +
  theme_bw(base_size = 18)       


