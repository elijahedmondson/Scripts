
GRSDF <- subset(GRSD, sex=='female')
GRSDM <- subset(GRSD, sex=='male')

GRSDF.IR <- subset(GRSD.IR, sex=='female')
GRSDM.IR <- subset(GRSD.IR, sex=='male')

names(GRSD)

attach(GRSD)

TAB = table(HCC, GROUP)
TAB

barplot(TAB, width = 1, xlab = "Irradiated", ylab = "Number of mice", beside = T, legend = c("No Disease", "Disease"), 
        args.legend = list(title = "x", x = "topright", cex = .9), ylim = c(0, 1200))



chisq.test(TAB, correct=T) #correct=T provides Yate's continuity correction
fisher.test(TAB, conf.int=T, conf.lev=0.95)

TAB



-barplot2(TAB, width = 1, xlab = "Irradiated", ylab = "Number of mice", beside = T, legend = c("No Disease", "Disease"), ylim = c(0, 1200))



