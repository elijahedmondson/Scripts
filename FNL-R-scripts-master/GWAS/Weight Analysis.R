library(HZE)
library(ggplot2)
library(doParallel)
library(foreach)
library(Rsamtools)
library(VariantAnnotation)
library(GenomicRanges)
library(regress)
library(MASS)
library(lmtest)
library(car)
library(DOQTL)

outdir = "~/Desktop/"
setwd(outdir)

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


addcovar = matrix(phenow$sex, ncol = 1, dimnames = list(rownames(phenow), "sex"))
cohort1 = subset(pheno, cohort == 1)
cohort2 = subset(pheno, cohort == 2)

HZE = subset(pheno, group == "HZE" & hardtumor == "1")
Gamma = subset(pheno, group == "Gamma" & hardtumor == "1")
Unirradiated = subset(pheno, group == "Unirradiated" & hardtumor == "1")

Male = subset(pheno, sex == 1)
Malew = Male[complete.cases(Male[,6]),]
Madd = matrix(Malew$sex, ncol = 1, dimnames = list(rownames(Malew), "sex"))

Fem = subset(pheno, sex == 0)
Femw = Fem[complete.cases(Fem[,6]),]
Fadd = matrix(Femw$sex, ncol = 1, dimnames = list(rownames(Femw), "sex"))



Weight.QTL = scanone.assoc(pheno = phenow, pheno.col = "weight", probs = model.probs, K = K, addcovar = addcovar, markers = MM_snps, sdp.file = sdp.file, ncl = 4)
Male.weight.QTL = scanone.assoc(pheno = Malew, pheno.col = "weight", probs = model.probs, K = K, addcovar = Madd, markers = MM_snps, sdp.file = sdp.file, ncl = 4)
Female.weight.QTL = scanone.assoc(pheno = Femw, pheno.col = "weight", probs = model.probs, K = K, addcovar = Fadd, markers = MM_snps, sdp.file = sdp.file, ncl = 4)

save(Weight.QTL, Male.weight.QTL, Female.weight.QTL, file = "~/Desktop/Weight.Corrected.QTL.MF.Rdata")

DOQTL:::plot.scanone.assoc(Weight.QTL, bin.size = 100, main = "Weight")
DOQTL:::plot.scanone.assoc(Male.weight.QTL, bin.size = 100, main = "Male.weight.QTL")
DOQTL:::plot.scanone.assoc(Female.weight.QTL, bin.size = 100, main = "Female.weight.QTL")

layout(matrix(2:1, nrow = 2, ncol = 1))
par(mfrow = c(1,1), mar=c(1, 4, 1, 1) + 0.1)
DOQTL:::plot.scanone.assoc(Male.weight.QTL, bin.size = 100, main = "Male.weight.QTL")
DOQTL:::plot.scanone.assoc(Female.weight.QTL, bin.size = 100, main = "Female.weight.QTL")

layout(matrix(1:6, nrow = 3, ncol = 2))
par(mfrow = c(3,2), mar=c(1, 4, 1, 1) + 0.1)
DOQTL:::plot.scanone.assoc(Weight.QTL, chr=1, bin.size = 50, main = "All Chr1", ylim=c(0,35))
DOQTL:::plot.scanone.assoc(Weight.QTL, chr=4, bin.size = 50, main = "All Chr4", ylim=c(0,35))
DOQTL:::plot.scanone.assoc(Male.weight.QTL, chr=1, bin.size = 50, main = "Male Chr1", ylim=c(0,35))
DOQTL:::plot.scanone.assoc(Male.weight.QTL, chr=4, bin.size = 50, main = "Male Chr4", ylim=c(0,35))
DOQTL:::plot.scanone.assoc(Female.weight.QTL, chr=1, bin.size = 50, main = "Female Chr1", ylim=c(0,35))
DOQTL:::plot.scanone.assoc(Female.weight.QTL, chr=4, bin.size = 50, main = "Female Chr4", ylim=c(0,35))


phenow <- pheno[complete.cases(pheno[,6]),]
phenow <- phenow[order(mean(pheno$days)),]





ggplot(pheno, aes(x = family, y = weight, color = cohort)) +
        geom_smooth(method=lm, fullrange = T) +
        geom_point(aes(color=factor(sex), shape=factor(sex)), size = 3) +
        scale_color_manual(values = c("red", "blue"), name="", labels = c("")) +
        scale_shape_manual(values = c(2,0), name="", labels = c("")) +
        theme_bw(base_size = 18) +
        ggtitle("") +
        theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"),
              legend.position=c(.9,.1))






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

p3 <- ggplot(Unirradiated, aes(x = reorder(family, days, FUN = median), y = days)) + geom_boxplot(notch = T, aes(fill = factor(family))) + geom_jitter() +
        theme_bw(base_size = 18) +
        ggtitle("Harderian Tumor Latency: Unirradiated") +
        xlab("HS/npt Family") +
        theme(axis.text = element_text(size = 14),
              legend.position = "none",
              panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_blank())
multiplot(p1,p2, cols = 1)





### Reorder
p1 <- ggplot(Unirradiated, aes(x = family, y = days)) + geom_boxplot(notch = T, aes(fill = factor(family))) + geom_jitter() +
        theme_bw(base_size = 18) +
        ggtitle("Overall Female Survival Time of HS/npt Mice by Family") +
        theme(legend.position = "none", axis.title.x = element_blank(), plot.margin=unit(c(1,1,1.5,1.2),"cm"))
p2 <- ggplot(Unirradiated, aes(x = reorder(family, days, FUN = median), y = days)) + geom_boxplot(notch = T, aes(fill = factor(family))) + geom_jitter() +
        theme_bw(base_size = 18) +
        ggtitle("Sorted Female Survival of HS/npt Mice by Family") +
        theme(legend.position = "none", axis.title.x = element_blank(), plot.margin=unit(c(1,1,1.5,1.2),"cm"))
multiplot(p1,p2, cols = 1)

p3 <- ggplot(cohort1w, aes(x = reorder(family, weight, FUN = median), y = weight)) +
        geom_boxplot(notch = T, aes(fill = factor(family))) + geom_jitter() +
        theme_bw(base_size = 18) +
        ggtitle("Cohort 1 Weight by Family") +
        theme(legend.position = "none", axis.title.x = element_blank(), plot.margin=unit(c(1,1,1.5,1.2),"cm"))

p4 <- ggplot(cohort2w, aes(x = reorder(family, weight, FUN = median), y = weight)) +
        geom_boxplot(notch = T, aes(fill = factor(family))) + geom_jitter() +
        theme_bw(base_size = 18) + ggtitle("Cohort 2 Weight by Family") +
        theme(legend.position = "none", axis.title.x = element_blank(), plot.margin=unit(c(1,1,1.5,1.2),"cm"))
multiplot(p3,p4, cols = 1)


p1 <- ggplot(pheno, aes(x = weight, fill = cohort)) +
        geom_histogram(colour="black", binwidth = 1) +
        facet_grid(cohort ~ .) +
        stat_function(fun=dbeta,args=fitdistr(phenow$weight,"beta",start=list(shape1=1,shape2=1))$estimate) +
        theme_bw() +
        ggtitle("Weight of HS/npt Mice by Family") +
        theme(legend.position = "none", axis.title.x = element_blank(), plot.margin=unit(c(1,1,1.5,1.2),"cm"))








#####HISTOGRAM#####
layout(matrix(3:1, 3, 1))
hist(pheno$days, breaks=150, col="blue", main = "")
lines(density(un$AML.Locus), col="black")
hist(pheno$days, breaks=150, col="green", main = "Gamma", xlab="Chromosome 2", prob = T, xlim = c(0, 182113224))
lines(density(gamma$AML.Locus), col="black")
hist(hze$AML.Locus, breaks=150, col="red", main = "HZE", xlab="Chromosome 2", prob = T, xlim = c(0, 182113224))
lines(density(hze$AML.Locus), col="black", lwd = 1)



AML.boot <- factor(pheno$days, levels = c("All.irradiated", "gamma", "HZE", "Unirradiated"),
                   labels = c("All Irradiated", "Gamma", "HZE", "Unirradiated"))
par(mfrow=c(1,1))
sm.density.compare(Thy.boot2$AML.Locus, Thy.boot2$TX, xlab = "Chromosome 2", lwd = 2.5)
title(main = "AML Adenoma: Resample Model Averaging")
colfill = c(2:(2+length(levels(AML.boot))))
legend("topright", levels(AML.boot), fill = colfill)
#sm.binomial.bootstrap()



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
        library(grid)

        # Make a list from the ... arguments and plotlist
        plots <- c(list(...), plotlist)

        numPlots = length(plots)

        # If layout is NULL, then use 'cols' to determine layout
        if (is.null(layout)) {
                # Make the panel
                # ncol: Number of columns of plots
                # nrow: Number of rows needed, calculated from # of cols
                layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                                 ncol = cols, nrow = ceiling(numPlots/cols))
        }

        if (numPlots==1) {
                print(plots[[1]])

        } else {
                # Set up the page
                grid.newpage()
                pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

                # Make each plot, in the correct location
                for (i in 1:numPlots) {
                        # Get the i,j matrix positions of the regions that contain this subplot
                        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

                        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                        layout.pos.col = matchidx$col))
                }
        }
}

