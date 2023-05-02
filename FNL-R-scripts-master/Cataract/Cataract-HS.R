###Cataract Analysis
library(readxl)
library(survival)
library(survminer)
library(mgcv)

data <- read_excel("C:/Users/edmondsonef/Desktop/Cataract/CATARACT_final.xlsx", sheet = "KM")#"CATARACT_fin")

fit <- survfit(Surv(data$Cat4.0day, data$Cat4.0) ~ data$groups, data = data)
# Visualize with survminer
ggsurvplot(fit, data = data, risk.table = TRUE, conf.int = T, pval= T)



Control <- subset(data, data$Exposure=='Control')
fit2 <- survfit(Surv(cat2.days, cat2) ~ Sex, data = Control)
ggsurvplot(fit2, data = Control, risk.table = TRUE, conf.int = T, pval= T)
print(fit2, print.rmean=TRUE)

rm(fit2)




fit <- survfit(Surv(days) ~ data$'groups', data = data)
# Visualize with survminer
ggsurvplot(fit, data = data, risk.table = TRUE, conf.int = F, pval= T)



newdata <- subset(data, data$`Parathyroid Adenoma`=='1')
fit2 <- survfit(Surv(days) ~ groups, data = newdata)
ggsurvplot(fit2, data = newdata, risk.table = TRUE, conf.int = F, pval= T)
print(fit2, print.rmean=TRUE)

rm(fit2)


### Univariate Cox regression ###
### Univariate Cox regression ###
### Univariate Cox regression ###
### Univariate Cox regression ###

res.cox <- coxph(Surv(data$cat2.days, data$cat2) ~ data$Exposure, data = data)
res.cox

covariates <- c(Sex, Exposure, Family, HZE.Ion)
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(data$cat2.days, data$cat2) ~ data$', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                       "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
HR <- as.data.frame(res.cox)

write.csv(res.cox, "C:/Users/edmondsonef/Desktop/catHR.csv")

### Multivariate Cox regression ###
### Multivariate Cox regression ###
### Multivariate Cox regression ###
### Multivariate Cox regression ###
data <- read_excel("C:/Users/edmondsonef/Desktop/Manuscripts/Cataract/GRSD.cataract.xlsx", sheet = "Harderian censored")
data <- read_excel("C:/Users/edmondsonef/Desktop/Manuscripts/Cataract/GRSD.cataract.xlsx", sheet = "NA censored")

model <- print.coxph(Surv(data$cat2.days, data$cat2) ~ Sex + Group, data =  data)
summary(model)

cox  <- model

# Prepare the columns
beta <- coef(cox)
se   <- sqrt(diag(cox$var))
p    <- 1 - pchisq((beta/se)^2, 1)
CI   <- round(confint(cox), 3)

# Bind columns together, and select desired rows
res <- cbind(beta, se = exp(beta), CI, p)
res <- res[c("age", "log(protime)"),]

write.csv(res, "C:/Users/edmondsonef/Desktop/catHR.csv")

ggforest(model, data = data, fontsize = 0.9)

###Testing the assumption of proportionality###
cox.zph(res.cox)
plot(cox.zph(res.cox))




### Multivariate Cox regression & Forrest Plot ###
### Multivariate Cox regression & Forrest Plot ###
### Multivariate Cox regression & Forrest Plot ###
### Multivariate Cox regression & Forrest Plot ###


ggforest(model, data = data, main = "Hazard Ratio: Merriam-Focht Score 2.0",
         cpositions = c(0.02, 0.15, 0.3), fontsize = 0.9,
         refLabel = "reference", noDigits = 2)


