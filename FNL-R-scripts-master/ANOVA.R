library(ggplot2)
library(ggpubr)
library(tidyverse)
library(broom)
library(AICcmodavg)

summary(data)

##Determine Normality
hist(data$Tubule)
hist(data$Inflammation)


##One-way ANOVA
one.way <- aov(Inflammation ~ Strain, data = data)
summary(one.way)


##Two-way ANOVA
two.way <- aov(Status ~ Tubule + Age, data = data)
summary(two.way)


##Adding interactions between variables
interaction <- aov(Inflammation ~ Age*Status, data = data)
summary(interaction)

##Adding blocking variable
blocking <- aov(Inflammation ~ Age + Status + block, data = data)
summary(blocking)


##Find the best-fit model
library(AICcmodavg)

model.set <- list(one.way, two.way, interaction, blocking)
model.names <- c("one.way", "two.way", "interaction", "blocking")

aictab(model.set, modnames = model.names)


##Check for homoscedasticity
par(mfrow=c(2,2))
plot(two.way)
par(mfrow=c(1,1))


##Do a post-hoc test
tukey.two.way<-TukeyHSD(two.way)
tukey.two.way


##Plot the results in a graph
tukey.plot.aov<-aov(Tubule ~ Inflammation:Status, data=data)

tukey.plot.test<-TukeyHSD(tukey.plot.aov)
plot(tukey.plot.test, las = 1)



#########Example: Reporting the results of ANOVA
###
###We found a statistically-significant difference in 
###average crop yield by both fertilizer type (f(2)=9.018, 
###p < 0.001) and by planting density (f(1)=15.316, p<0.001).
###
###A Tukey post-hoc test revealed that fertilizer mix 3 
###resulted in a higher yield on average than fertilizer mix 1 
###(0.59 bushels/acre), and a higher yield on average than 
###fertilizer mix 2 (0.42 bushels/acre). Planting density was 
###also significant, with planting density 2 resulting in an 
###higher yield on average of 0.46 bushels/acre over planting 
###density 1.
###
###A subsequent groupwise comparison showed the strongest yield 
###gains at planting density 2, fertilizer mix 3, suggesting that 
###this mix of treatments was most advantageous for crop growth 
###under our experimental conditions.




