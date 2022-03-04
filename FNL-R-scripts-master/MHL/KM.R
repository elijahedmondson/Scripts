#data <- read_excel("JBM 336.xlsx")

library(survMisc)
library(survival)
library(jskm)
library(survey)
library(survminer)
library(survival)
library(ggplot2)
library(ggfortify)

####### OPTION 1 ####### 
####### OPTION 1 ####### 
####### OPTION 1 ####### 
####### OPTION 1 ####### 
####### OPTION 1 ####### 
####### OPTION 1 ####### 
####### OPTION 1 ####### 
####### OPTION 1 ####### 
####### OPTION 1 ####### 

fit <- survfit(Surv(Days)~Group, data=data)

jskm(fit)
jskm(fit, ci = F, cumhaz = F, legendposition = c(0.2,0.2),  mark = F, table = T, ylab = "Cumulative incidence (%)", surv.scale = "percent")#, pval =T, pval.size = 6, pval.coord = c(300, 0.7))

tiff("JBM336.tiff", units="in", width=10, height=5, res=300)
jskm(fit, ci = F, cumhaz = F, legendposition = c(0.2,0.2),  mark = F, ylab = "Cumulative incidence (%)", surv.scale = "percent")
dev.off()


####### OPTION 2 ####### 
####### OPTION 2 ####### 
####### OPTION 2 ####### 
####### OPTION 2 ####### 
####### OPTION 2 ####### 
####### OPTION 2 #######
####### OPTION 2 ####### 
####### OPTION 2 ####### 
####### OPTION 2 #######


model_fit <- survfit(Surv(Age, Status) ~ Groups, data = data)

autoplot(model_fit) + 
  labs(x = "\n Survival Time (Days) ", y = "Survival Probabilities \n", 
       title = "Survival Times \n JBM 336 \n") + 
  theme_bw()

autoplot(survfit(Surv(Age, Status)~ data$'CNS Polyglucosan Body', data = data), conf.int = FALSE, censor = T)
autoplot(survfit(Surv(Age, Status) ~ Groups, data = data))
autoplot(survfit(Surv(Age, Status) ~ Groups, data = data), conf.int = FALSE, censor = T)


####### BOX COX Transform ####### 
####### BOX COX Transform ####### 
####### BOX COX Transform ####### 
####### BOX COX Transform ####### 
####### BOX COX Transform ####### 
####### BOX COX Transform ####### 

library(MASS)

#fit linear regression model
model <- lm(data$Age~data$Status)
plot(model)

bc <- boxcox(data$Age~data$Status)
(lambda <- bc$x[which.max(bc$y)])

new_model <- lm(((data$Age^lambda-1)/lambda) ~ data$Status)
