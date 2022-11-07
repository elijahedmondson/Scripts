library(car)
library(HZE)
library(ggplot2)

GRSD.pheno <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = GRSD.pheno$row.names, sex = as.numeric(GRSD.pheno$sex == "M"),
                   cohort = as.numeric(GRSD.pheno$Cohort),
                   days = as.numeric(GRSD.pheno$days),
                   catdays = as.numeric(GRSD.pheno$Cataract.2.0.Score.Days),
                   cataract = as.numeric(GRSD.pheno$Cataract.2.0.Score.Event),
                   BCS = as.character(GRSD.pheno$BCS),
                   group = as.character(GRSD.pheno$groups),
                   family = as.numeric(GRSD.pheno$family),
                   weight.old = as.numeric(GRSD.pheno$Mass),
                   weight = as.numeric(GRSD.pheno$Weight.corrected),
                   hardtumor = as.numeric(GRSD.pheno$Harderian.Tumor))
All = subset(pheno, hardtumor == "1")
HZE = subset(pheno, group == "HZE" & hardtumor == "1")
Gamma = subset(pheno, group == "Gamma" & hardtumor == "1")
Unirradiated = subset(pheno, group == "Unirradiated" & hardtumor == "1")


#Harderian Tumor by Family and Days
p1 <- ggplot(HZE, aes(x = reorder(family, days, FUN = median), y = days)) + geom_boxplot(notch = T, aes(fill = factor(family))) + geom_jitter() +
        theme_bw(base_size = 18) +
        ggtitle("Harderian Tumor Latency: HZE Ion Irradiated") +
        xlab("HS/npt Family") +
        theme(axis.text = element_text(size = 14),
              legend.position = "none",
              panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_blank())

p2 <- ggplot(Gamma, aes(x = reorder(family, days, FUN = median), y = days)) + geom_boxplot(notch = T, aes(fill = factor(family))) + geom_jitter() +
        theme_bw(base_size = 18) +
        ggtitle("Harderian Tumor Latency: Gamma-ray Irradiated") +
        xlab("HS/npt Family") +
        theme(axis.text = element_text(size = 14),
              legend.position = "none",
              panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_blank())

p3 <- ggplot(All, aes(x = reorder(family, days, FUN = median), y = days)) + geom_boxplot(notch = T, aes(fill = factor(family))) + geom_jitter() +
        theme_bw(base_size = 18) +
        ggtitle("Harderian Tumor Latency: Unirradiated") +
        xlab("HS/npt Family") +
        theme(axis.text = element_text(size = 14),
              legend.position = "none",
              panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_blank())
multiplot(p1,p2, cols = 1)

ggplot(All, aes(x = reorder(group, days, FUN = median), y = days)) + geom_boxplot(notch = T, aes(fill = factor(group))) + geom_jitter() +
        theme_bw(base_size = 18) +
        ggtitle("Harderian Tumor Latency: Unirradiated") +
        xlab("HS/npt Family") +
        theme(axis.text = element_text(size = 14),
              legend.position = "none",
              panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_blank())


##### Type II SS #####
model <- anova(lm(data = pheno, cataract ~ family * group * cohort))
ss <- model$"Sum Sq"
model <- cbind(model, PctExp = ss/sum(ss)*100)
model

##### Type III SS #####
options(contrasts = c("contr.sum","contr.poly"))
model <- lm(data = All, days ~ group * family * cohort)
drop1(model, .~., test="F")

af <- anova(fit)
afss <- af$"Sum Sq"
print(cbind(af,PctExp=afss/sum(afss)*100))
