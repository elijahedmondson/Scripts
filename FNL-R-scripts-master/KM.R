
###Everything after a "#" is a developer comment
###
###The first time you work with R, you will need to download all your packages

install.packages("survival")

###
###Do this for all packages... the second time you open R, you can start with the following
###

###
###Start by loading your packages
###
library(plyr)
library(dplyr)
library(readr)
library(survival)
library(survminer)
library(readxl)
library(ggplot2)
library(gridExtra)
library(readxl)
library(ggpubr)
library(Rmisc)

###
### Then load in your file
###

###
###Then create your survival model
###
surv_object = Surv(time = data$Days, event = data$Group)

# Regress on a constant
fit <- survfit(surv_object ~ 1)

# Plot the fit
ggsurvplot(fit, data.frame(time = data$Days, event = data$Group), conf.int=FALSE)

sarcoma.fit <- survfit(Surv(data$Days, data$Sex) ~ pheno$'Group', data = data)
HN.fit <- survfit(Surv(data$Days, pheno$Hematopoietic) ~ pheno$'Group', data = pheno)

###
### Visualize with survminer
###
ggsurvplot(sarcoma.fit, data = pheno, risk.table = TRUE, conf.int = F, pval= T)
ggsurvplot(HN.fit, data = pheno, risk.table = TRUE, conf.int = F, pval= T)







###
### Other code
###
###
### Other code
###
###
### Other code
###
###
### Other code
###


fit2 <- survfit(Surv(Days, data$'Pulmonary carcinoma') ~ Groups, data = data)
ggsurvplot(fit2, data = data, risk.table = TRUE, conf.int = F, pval= T)
print(fit2, print.rmean=TRUE)

rm(fit2)



library(survival)
library(survminer)

fit <- survfit(Surv(Age) ~ data$'Group', data = data)
# Visualize with survminer
ggsurvplot(fit, data = data, risk.table = TRUE, conf.int = F, pval= T)



newdata <- subset(data, data$'Proliferative Pulmonary Lesions'=='1')
fit2 <- survfit(Surv(Days) ~ Groups, data = newdata)
ggsurvplot(fit2, data = newdata, risk.table = TRUE, conf.int = F, pval= T)
print(fit2, print.rmean=TRUE)

rm(fit2)
