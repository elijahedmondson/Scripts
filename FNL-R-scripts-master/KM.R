
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
pheno <- read_excel("C:/Users/edmondsonef/Desktop/16-306 Animal Data.xlsx")

###
###Then create your survival model
###

sarcoma.fit <- survfit(Surv(pheno$age, pheno$Sarcoma) ~ pheno$'Group', data = pheno)
HN.fit <- survfit(Surv(pheno$age, pheno$Hematopoietic) ~ pheno$'Group', data = pheno)

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

newdata <- subset(pheno, pheno$`Parathyroid Adenoma`=='1')
fit2 <- survfit(Surv(days) ~ groups, data = newdata)
ggsurvplot(fit2, data = newdata, risk.table = TRUE, conf.int = F, pval= T)
print(fit2, print.rmean=TRUE)

rm(fit2)



library(survival)
library(survminer)

fit <- survfit(Surv(Age) ~ data$'Group', data = data)
# Visualize with survminer
ggsurvplot(fit, data = data, risk.table = TRUE, conf.int = F, pval= T)



newdata <- subset(data, data$`Parathyroid Adenoma`=='1')
fit2 <- survfit(Surv(days) ~ groups, data = newdata)
ggsurvplot(fit2, data = newdata, risk.table = TRUE, conf.int = F, pval= T)
print(fit2, print.rmean=TRUE)

rm(fit2)
