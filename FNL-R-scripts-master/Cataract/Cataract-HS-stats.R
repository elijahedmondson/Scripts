library(readxl)
library(tidyverse)
library(DescTools)
library(lme4)
library(broom.mixed)
library(kableExtra)
library(xtable)
library(emmeans)
library(ggsci)
library(superdiag)
library(mcmcplots)
library(gridExtra)
library(plyr)
library(forcats)
library(gghighlight)
library(dplyr)
library(ggplot2)
library(readxl)
library(survival)
library(survminer)
library(flexsurv)
library(survtools)
library(finalfit)
library(gtsummary)
library(tidyverse)
library(ggfortify)
library(ggplot2)
library(ggspectra)
library(ggrepel)
library(survival)
library(ranger)
library(dplyr)
library(ggfortify)

library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)

#data <- read_excel(path = "C:/Users/edmondsonef/Desktop/Cataract/CAT_final.xlsx")
#data <- read_excel(path = "C:/Users/edmondsonef/Desktop/Cataract/CAT_final.xlsx", sheet = "RM_HardAd_DND")
data <- read_excel(path = "C:/Users/edmondsonef/Desktop/Cataract/CAT_final.xlsx", sheet = "RM_DND")
data <- data %>%  mutate(Family = as.character(family), 
                         BCS = as.ordered(BCS),
                         Sex = sex,
                         Treatment = relevel(as.factor(groups), ref = "Unirradiated"),
                         `Number of Tumors` = relevel(as.factor(`Number of Tumors`), ref = "No tumors"))
#data <- data %>%  filter(Age >= 552)

#####Plots by Family
#####
cats$families <-  as.character(cats$Family)
#highlight_df <- data %>% filter(data$family == c(6, 36, 50, 47, 20, 34))


mu <- ddply(cats, "Treatment", summarise, grp.mean=mean(Cat_score))
ggplot(cats, aes(x=Cat_score, color=Sex)) +
  geom_histogram(fill="white", alpha=0.5, position="identity", bins =200)+
  geom_vline(data=mu, aes(xintercept=grp.mean),
             linetype="dashed")+
  theme(legend.position="top")+
  facet_grid(Treatment ~ .)


cats %>%
  mutate(family = fct_reorder(families, Family)) %>%
  ggplot(aes(group = Family, x = Family, y = Cat_score, fill = families)) +
  geom_boxplot() +
  stat_summary(fun = "mean", geom = "point", shape = 8, size = 2, color = "white")+
  theme_bw() +
  theme(legend.position = "none")+
  facet_grid(Sex ~ .)



plot1 <- cats %>%
  mutate(Family = fct_reorder(families, Cat_score)) %>%
  ggplot(aes(group = Family, x = Family, y = Cat_score, fill = families)) +
  geom_boxplot() +
  #stat_summary(fun = "mean", geom = "point", shape = 8, size = 2, color = "white")+
  theme_bw() +
  theme(legend.position = "none")+
  facet_grid(Treatment ~ .)

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Cat_score_by_treatment.tiff", units="in", width=10, height=6, res=200)
plot1
dev.off()



plot2 <- cats %>%
  mutate(Family = fct_reorder(families, Cat_score)) %>%
  ggplot(aes(group = Family, x = Family, y = Cat_score, fill = families)) +
  geom_boxplot() +
  #stat_summary(fun = "mean", geom = "point", shape = 8, size = 2, color = "white")+
  theme_bw() +
  theme(legend.position = "none")+
  facet_grid(Sex ~ .)

setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Cat_score_by_sex.tiff", units="in", width=10, height=6, res=200)
plot2
dev.off()



cats %>%
  ggplot(aes(group = Treatment, x = Treatment, y = Cat_score, fill = Treatment)) +
  geom_boxplot() +
  stat_summary(fun = "mean", geom = "point", shape = 8, size = 2, color = "white")+
  theme_bw() +
  theme(legend.position = "none")+
  facet_grid(Sex ~ .)


#####

#####KM by Radiation
#####
Status = data$`Cat 0.5`
Day = data$`Cat 0.5 day`
model_fit <- survfit(Surv(Day, Status) ~ groups, data = data)
Cat0.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 0.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Groups",legend.labs=c("Gamma","HZE","Unirradiated"))
Status = data$`Cat 1.0`
Day = data$`Cat 1.0 day`
model_fit <- survfit(Surv(Day, Status) ~ groups, data = data)
Cat1.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 1.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Groups",legend.labs=c("Gamma","HZE","Unirradiated"))
Status = data$`Cat 1.5`
Day = data$`Cat 1.5 day`
model_fit <- survfit(Surv(Day, Status) ~ groups, data = data)
Cat1.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 1.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Groups",legend.labs=c("Gamma","HZE","Unirradiated"))
Status = data$`Cat 2.0`
Day = data$`Cat 2.0 day`
model_fit <- survfit(Surv(Day, Status) ~ groups, data = data)
Cat2.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 2.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Groups",legend.labs=c("Gamma","HZE","Unirradiated"))
Status = data$`Cat 2.5`
Day = data$`Cat 2.5 day`
model_fit <- survfit(Surv(Day, Status) ~ groups, data = data)
Cat2.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 2.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Groups",legend.labs=c("Gamma","HZE","Unirradiated"))
Status = data$`Cat 3.0`
Day = data$`Cat 3.0 day`
model_fit <- survfit(Surv(Day, Status) ~ groups, data = data)
Cat3.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 3.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Groups",legend.labs=c("Gamma","HZE","Unirradiated"))
Status = data$`Cat 3.5`
Day = data$`Cat 3.5 day`
model_fit <- survfit(Surv(Day, Status) ~ groups, data = data)
Cat3.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 3.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Groups",legend.labs=c("Gamma","HZE","Unirradiated"))
Status = data$`Cat 4.0`
Day = data$`Cat 4.0 day`
model_fit <- survfit(Surv(Day, Status) ~ groups, data = data)
Cat4.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 4.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Groups",legend.labs=c("Gamma","HZE","Unirradiated"))

# List of ggsurvplots
require("survminer")
splots <- list()
splots[[1]] <- Cat0.5
splots[[2]] <- Cat1.0
splots[[3]] <- Cat1.5
splots[[4]] <- Cat2.0
splots[[5]] <- Cat2.5
splots[[6]] <- Cat3.0
splots[[7]] <- Cat3.5
splots[[8]] <- Cat4.0


setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Plot.tiff", units="in", width=12, height=15, res=200)
arrange_ggsurvplots(splots, print = TRUE, ncol = 2, nrow = 4, risk.table.height = 0.4)
dev.off()
#####

#####KM by Sex
#####
Status = data$`Cat 0.5`
Day = data$`Cat 0.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat0.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 0.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 1.0`
Day = data$`Cat 1.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat1.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 1.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 1.5`
Day = data$`Cat 1.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat1.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 1.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 2.0`
Day = data$`Cat 2.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat2.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 2.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 2.5`
Day = data$`Cat 2.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat2.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 2.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 3.0`
Day = data$`Cat 3.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat3.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 3.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 3.5`
Day = data$`Cat 3.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat3.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 3.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 4.0`
Day = data$`Cat 4.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat4.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 4.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")

# List of ggsurvplots
require("survminer")
splots <- list()
splots[[1]] <- Cat0.5
splots[[2]] <- Cat1.0
splots[[3]] <- Cat1.5
splots[[4]] <- Cat2.0
splots[[5]] <- Cat2.5
splots[[6]] <- Cat3.0
splots[[7]] <- Cat3.5
splots[[8]] <- Cat4.0


setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Plot.tiff", units="in", width=12, height=15, res=200)
arrange_ggsurvplots(splots, print = TRUE, ncol = 2, nrow = 4, risk.table.height = 0.4)
dev.off()
#KM plots

#####KM by sex within radiation

#####KM Unirradiated by sex

#####KM Unirradiated by sex
#####
data <- data %>%  filter(Treatment == "Unirradiated")

Status = data$`Cat 0.5`
Day = data$`Cat 0.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat0.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 0.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 1.0`
Day = data$`Cat 1.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat1.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 1.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 1.5`
Day = data$`Cat 1.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat1.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 1.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 2.0`
Day = data$`Cat 2.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat2.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 2.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 2.5`
Day = data$`Cat 2.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat2.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 2.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 3.0`
Day = data$`Cat 3.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat3.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 3.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 3.5`
Day = data$`Cat 3.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat3.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 3.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 4.0`
Day = data$`Cat 4.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat4.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 4.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")

# List of ggsurvplots
require("survminer")
splots <- list()
splots[[1]] <- Cat0.5
splots[[2]] <- Cat1.0
splots[[3]] <- Cat1.5
splots[[4]] <- Cat2.0
splots[[5]] <- Cat2.5
splots[[6]] <- Cat3.0
splots[[7]] <- Cat3.5
splots[[8]] <- Cat4.0


setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Plot.tiff", units="in", width=12, height=15, res=200)
arrange_ggsurvplots(splots, print = TRUE, ncol = 2, nrow = 4, risk.table.height = 0.4)
dev.off()


#####

#####
#####KM HZE by sex
#####
data <- data %>%  filter(Treatment == "HZE")

Status = data$`Cat 0.5`
Day = data$`Cat 0.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat0.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 0.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 1.0`
Day = data$`Cat 1.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat1.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 1.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 1.5`
Day = data$`Cat 1.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat1.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 1.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 2.0`
Day = data$`Cat 2.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat2.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 2.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 2.5`
Day = data$`Cat 2.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat2.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 2.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 3.0`
Day = data$`Cat 3.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat3.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 3.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 3.5`
Day = data$`Cat 3.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat3.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 3.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 4.0`
Day = data$`Cat 4.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat4.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 4.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")

# List of ggsurvplots
require("survminer")
splots <- list()
splots[[1]] <- Cat0.5
splots[[2]] <- Cat1.0
splots[[3]] <- Cat1.5
splots[[4]] <- Cat2.0
splots[[5]] <- Cat2.5
splots[[6]] <- Cat3.0
splots[[7]] <- Cat3.5
splots[[8]] <- Cat4.0


setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Plot.tiff", units="in", width=12, height=15, res=200)
arrange_ggsurvplots(splots, print = TRUE, ncol = 2, nrow = 4, risk.table.height = 0.4)
dev.off()

#####
#####KM Gamma by sex
#####
data <- data %>%  filter(Treatment == "Gamma")

Status = data$`Cat 0.5`
Day = data$`Cat 0.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat0.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 0.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 1.0`
Day = data$`Cat 1.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat1.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 1.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 1.5`
Day = data$`Cat 1.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat1.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 1.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 2.0`
Day = data$`Cat 2.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat2.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 2.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 2.5`
Day = data$`Cat 2.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat2.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 2.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 3.0`
Day = data$`Cat 3.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat3.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 3.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 3.5`
Day = data$`Cat 3.5 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat3.5 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 3.5",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")
Status = data$`Cat 4.0`
Day = data$`Cat 4.0 day`
model_fit <- survfit(Surv(Day, Status) ~ sex, data = data)
Cat4.0 <- ggsurvplot2(model_fit, censor = F, data=data, xlab = "Days", title="Score 4.0",
                      conf.int = T, pval = T, risk.table = F, 
                      legend="right",legend.title="Sex")

# List of ggsurvplots
require("survminer")
splots <- list()
splots[[1]] <- Cat0.5
splots[[2]] <- Cat1.0
splots[[3]] <- Cat1.5
splots[[4]] <- Cat2.0
splots[[5]] <- Cat2.5
splots[[6]] <- Cat3.0
splots[[7]] <- Cat3.5
splots[[8]] <- Cat4.0


setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("Plot.tiff", units="in", width=12, height=15, res=200)
arrange_ggsurvplots(splots, print = TRUE, ncol = 2, nrow = 4, risk.table.height = 0.4)
dev.off()
#####



##### 
##### CoxPH risk estimates
##### 

#data <- read_excel(path = "C:/Users/edmondsonef/Desktop/Cataract/CAT_final.xlsx")
#data <- read_excel(path = "C:/Users/edmondsonef/Desktop/Cataract/CAT_final.xlsx", sheet = "RM_HardAd_DND")
data <- read_excel(path = "C:/Users/edmondsonef/Desktop/Cataract/CAT_final.xlsx", sheet = "RM_DND")
data <- data %>%  mutate(Family = as.character(family), 
                         BCS = as.ordered(BCS),
                         Sex = sex,
                         Treatment = relevel(as.factor(groups), ref = "Unirradiated"),
                         `Number of Tumors` = relevel(as.factor(`Number of Tumors`), ref = "No tumors"))

Status = data$`Cat 2.0`
Day = data$`Cat 2.0 day`

coxfit <- coxph(Surv(Day, Status) ~ Treatment + Sex, data = data)
summary(coxfit)
ggforest(coxfit, data = data, main = "Hazard Ratio: Merriam-Focht Score 2.0",
         cpositions = c(0.02, 0.15, 0.3), fontsize = 1,
         refLabel = "reference", noDigits = 3)


coxfit <- coxph(Surv(Day, Status) ~ Treatment + Sex + Treatment:Sex, data = data)
summary(coxfit)

coxfit <- coxph(Surv(Day, Status) ~ Treatment + Sex + `Number of Tumors`, data = data)
summary(coxfit)
ggforest(coxfit, data = data, main = "Hazard Ratio: Merriam-Focht Score 2.0",
         cpositions = c(0.02, 0.15, 0.3), fontsize = 1,
         refLabel = "reference", noDigits = 3)



coxfit <- coxph(Surv(Day, Status) ~ `Number of Tumors`, data = data)
summary(coxfit)
ggforest(coxfit, data = data, main = "Hazard Ratio: Merriam-Focht Score 2.0",
         cpositions = c(0.02, 0.15, 0.3), fontsize = 1,
         refLabel = "reference", noDigits = 3)


# aa_fit <-aareg(Surv(Day, Status) ~ Family, data = data)
# autoplot(aa_fit, ncol = 2)+
#   theme_bw() +
#   theme(legend.position = "none")

##### 

#####  Test effects
##### 
#Use a likelihood ratio to test interaction of radiation and sex

Status = data$`Cat 2.0`
Day = data$`Cat 2.0 day`
coxfit1 <- coxph(Surv(Day, Status) ~ Treatment + Sex, data = data)
coxfit2 <- coxph(Surv(Day, Status) ~ Treatment + Sex + Treatment:Sex, data = data)
anova(coxfit1,coxfit2,test="Chisq")


#####


#####
#data <- read_excel(path = "C:/Users/edmondsonef/Desktop/Cataract/CAT_final.xlsx")
#data <- read_excel(path = "C:/Users/edmondsonef/Desktop/Cataract/CAT_final.xlsx", sheet = "RM_HardAd_DND")
data <- read_excel(path = "C:/Users/edmondsonef/Desktop/Cataract/CAT_final.xlsx", sheet = "RM_DND")
data <- data %>%  mutate(Family = as.character(family), 
                         BCS = as.ordered(BCS),
                         Sex = sex,
                         Treatment = relevel(as.factor(groups), ref = "Unirradiated"),
                         `Number of Tumors` = relevel(as.factor(`Number of Tumors`), ref = "No tumors"))

Status = data$`Cat 2.0`
Day = data$`Cat 2.0 day`

covariate_names <- c(Treatment = relevel(as.factor(data$groups), ref = "Unirradiated"),
                     sex=data$sex,
                     `sex:female`="F",
                     Family = as.numeric(data$family),
                     CoatColor=data$"coat color",
                     BCS = as.ordered(data$BCS),
                     Amyloidosis=data$Amyloidosis,
                     AML=data$AML,
                     GCT=data$GCT,
                     Hard_ACA=data$Hard_ACA,
                     Hard_Ad=data$Hard_Ad,
                     HCC=data$HCC,
                     Hist_Sarc=data$Hist_Sarc,
                     LSA_PreT=data$LSA_PreT,
                     `LSA_B-cell`=data$`LSA_B-cell`,
                     Mamm_ACA=data$Mamm_ACA,
                     OSA=data$OSA,
                     Pulmonary_ACA = data$`Pulmonary Adenocarcinoma`,
                     Pituitary_Ad=data$Pituitary_Ad,
                     STS=data$STS,
                     Thyroid =data$Thyroid,
                     Pigment_Dispersion = data$Pigment_Dispersion,
                     Did_not_dilate = data$`Did Not Dilate`)


#Univariate Analyses
# It is common practice to perform univariate analyses of all covariates first and take only those into the multivariate analysis which were significant to some level in the univariate analysis. (I see some pros and strong cons with this procedure, but am open to learn more on this). The univariate part can easily be achieved using purrr's map function. A forest plot, as already said, will happily plot multiple results, even if they come as a list.

df <- data %>% mutate(sex=rename_factor(sex, `M` = "male", `F` = "female"))

map(vars(Treatment, sex, #Family,
         #`coat color`,
         Amyloidosis, Pigment_Dispersion, `Did Not Dilate`, AML, GCT,
         Hard_ACA, Hard_Ad, HCC, Hist_Sarc, LSA_PreT,
         `LSA_B-cell`,Mamm_ACA,OSA,`Pulmonary Adenocarcinoma`,
         Pituitary_Ad,STS,Thyroid,
         `Number of Tumors`), function(by)
{
  analyse_multivariate(df,
                       vars(Day, Status),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="OS"),
              orderer = ~order(HR),
              labels_displayed = c("endpoint", "factor", "n"),
              ggtheme = ggplot2::theme_bw(base_size = 12))


map(vars(Family), function(by)
         {
           analyse_multivariate(df,
                                vars(Day, Status),
                                covariates = list(by), # covariates expects a list
                                covariate_name_dict = covariate_names)
         }) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="OS"),
              orderer = ~order(HR),
              labels_displayed = c("endpoint", "factor", "n"),
              ggtheme = ggplot2::theme_bw(base_size = 12))


###Multivariate analysis

data %>%
  mutate(sex=rename_factor(sex, `M` = "male", `F` = "female")) %>%
  analyse_multivariate(vars(Day, Status),
                       covariates = vars(Treatment, `Number of Tumors`), 
                       covariate_name_dict = covariate_names) ->
  result
print(result)


forest_plot(result,
            factor_labeller = covariate_names,
            endpoint_labeller = c(time="OS"),
            orderer = ~order(HR),
            labels_displayed = c("endpoint", "factor", "n"),
            ggtheme = ggplot2::theme_bw(base_size = 10),
            relative_widths = c(1, 1.5, 1),
            HR_x_breaks = c(0.25, 0.5, 0.75, 1, 1.5, 2))

