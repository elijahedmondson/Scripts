# install.packages("survival")
library(splines)
library(survival)
library(KMsurv)
library(OIsurv)
library(rms)
library(ggplot2)
library(GGally)


require(ggplot2)
require(survival)
HCC.surv <- survfit(Surv(pheno1$days, pheno1$HCC) ~ pheno1$HCC.translocation)
ggsurv(HCC.surv)


Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Total-Table 1.csv")
pheno = data.frame(row.names = Total$corrected, original = Total$original,
                   group = Total$group,
                   sex = as.numeric(Total$sex == "M"),
                   days = as.numeric(Total$days),
                   OSA = as.numeric(Total$Osteosarcoma),
                   Mets = as.numeric(Total$Metastatic.Tumors))

Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = Total$row.names, rownames = Total$row.names,
                   sex = as.numeric(Total$sex == "M"),
                   days = as.numeric(Total$days),
                   HCC...translocation = as.numeric(Total$HCC...translocation),
                   HCC.translocation = as.numeric(Total$HCC.translocation),
                   HCC = as.numeric(Total$Hepatocellular.Carcinoma))

pheno1 = pheno[which(Total$Hepatocellular.Carcinoma=="1"),]
pheno1 = na.omit(pheno1)


head(pheno1)


time <- pheno1$days
event <- pheno1$HCC
pheno1$SurvObj <- Surv(time, event)

km.as.one <- npsurv(pheno1$SurvObj ~ 1, conf.type = "log-log")
survplot(km.as.one, title = "All Cases")
km.by.sex <- npsurv(pheno1$SurvObj ~ pheno1$sex, conf.type = "log-log")
survplot(km.by.sex)
km.by.pcr <- npsurv(pheno1$SurvObj ~ pheno1$HCC.translocation, conf.type = "log-log")
survplot(km.by.pcr)
km.by.pcr.amount <- npsurv(pheno1$SurvObj ~ pheno1$HCC...translocation, conf.type = "log-log")
ggsurv(km.by.pcr.amount)

survplot(Surv(pheno2$days, pheno2$HCC) ~ pheno2$HCC...translocation, data = pheno2)

OSA <- pheno[which(pheno$OSA == 1), ]
head(OSA)
OSA
allele <- c(1,0,1,2,2,1,0,1,1,1,1,2,1,1,1,1,2,1,2)

OSA <- cbind(OSA, allele)
time <- OSA$days
event <- (OSA$OSA)
OSA$SurvObj <- with(OSA, Surv(time, event))
OSA

km.as.one <- npsurv(OSA$SurvObj ~ 1, conf.type = "log-log")
km.by.sex <- npsurv(OSA$SurvObj ~ OSA$sex, conf.type = "log-log")
km.by.mets <- npsurv(OSA$SurvObj ~ OSA$Mets, conf.type = "log-log")
km.by.allele <- npsurv(OSA$SurvObj ~ OSA$allele, conf.type = "log-log")

survplot(km.as.one, title = "All Osteosarcoma Cases")

summary(km.by.allele)
survplot(km.by.allele, stitle = 'Osteosarcoma by Allele', xlab = 'Time (days)', 
         label.curves = F, conf = "bands", levels.only  = T, 
         n.risk = F, y.n.risk = -0, cex.n.risk = 0.7)
legend(40, .85, c("No BALB/cJ Allele", "BALB/cJ Heterozygous", "BALB/cJ Homozygous") , lty=c(1,3) )


survplot(km.by.mets, stitle = 'Osteosarcoma by Allele', xlab = 'Time (days)', n.risk = T)
survplot(km.by.sex)

#Log-logistic parametric model coefficients
loglogistic <- survreg(Surv(time, event) ~ OSA, dist="loglogistic")
summary(loglogistic)

# Cox proportional hazard model - coefficients and hazard rates
coxph <- coxph(Surv(time,event) ~ X, method="breslow")
summary(coxph)








# Descriptive statistics
summary(time)
summary(event)
summary(X)
summary(group)

# Kaplan-Meier non-parametric analysis
kmsurvival <- survfit(Surv(time,event) ~ 1)
summary(kmsurvival)
plot(kmsurvival, mark.time=TRUE, mark=1, col=1, lty=1, lwd=3, cex=1, log=FALSE, xlab="Days Post-Irradiation", ylab="Survival Probability")
title("GRSD Lymphoma")

# Kaplan-Meier non-parametric analysis by group
kmsurvival1 <- survfit(Surv(time, event) ~ group.broad)
summary(kmsurvival1)
plot(kmsurvival1, conf.int="both", mark.time=TRUE, 
     mark=4, col=1, lty=1:5, lwd=2, cex=1, log=FALSE, xscale=1, yscale=1,  
     firstx=0, firsty=1, ymin=0, xlab="Days Post-Irradiation", ylab="Survival Probability")
legend(30, .3, c("Gamma", "HZE", "Unirradiated"), lty = 1:5, lwd=2) 
title("GRSD Lymphoma")

#The previous without 95% CI
kmsurvival1 <- survfit(Surv(time, event) ~ group.broad)
summary(kmsurvival1)
plot(kmsurvival1, conf.int="none", mark.time=TRUE, 
     mark=4, col=1, lty=1:5, lwd=2, cex=1, log=FALSE, xscale=1, yscale=1,  
     firstx=0, firsty=1, ymin=0, xlab="Days Post-Irradiation", ylab="Survival Probability")
legend(30, .3, c("Gamma", "HZE", "Unirradiated"), lty = 1:5, lwd=2) 
title("GRSD Lymphoma")

# Nelson-Aalen non-parametric analysis
nasurvival <- survfit(coxph(Surv(time,event)~1), type="aalen")
summary(nasurvival)
plot(nasurvival, xlab="Time", ylab="Survival Probability")


# Cox proportional hazard model - coefficients and hazard rates
coxph <- coxph(Surv(time,event) ~ X, method="breslow")
summary(coxph)


# Exponential, Weibull, and log-logistic parametric model coefficients
exponential <- survreg(Surv(time,event) ~ X, dist="exponential")
summary(exponential)

weibull <- survreg(Surv(time,event) ~ X, dist="weibull")
summary(weibull)

loglogistic <- survreg(Surv(time,event) ~ X, dist="loglogistic")
summary(loglogistic)

