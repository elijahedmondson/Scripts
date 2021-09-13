library(survival)
library(haven)
Y= Surv(time = data$Age, event = data$Event)
kmfit = survfit(Y ~ data$Group)
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
surv_object = Surv(time = data$Age, event = data$Event)

# Regress on a constant
fit <- survfit(surv_object ~ 1)

# Plot the fit
ggsurvplot(fit, data.frame(time = data$Age, event = data$Group), conf.int=FALSE)

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
