### Nonparametric Bootstrap Resampling with Replacement ###

library(BSgenome.Mmusculus.UCSC.mm10)
library(doParallel)
library(VariantAnnotation)
library(survival)
library(GenomicRanges)
library(regress)
library(MASS)
library(DOQTL)
library(lmtest)
library(HZE)
library(dplyr)
library(sm)
options(stringsAsFactors = F)
load(file = "~/Desktop/R/QTL/WD/GRSD.Rdata")
load("/Users/elijah/Desktop/R/QTL/WD/hs.colors.Rdata")
setwd("~/Desktop/files")
outdir = "~/Desktop/files"
Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = Total$row.names, rownames = Total$row.names,
                   family = as.numeric(Total$family),
                   sex = as.numeric(Total$sex == "M"),
                   cohort = as.numeric(Total$Cohort),
                   group = as.character(Total$groups),
                   unirradiated = as.numeric(Total$Unirradiated),
                   days = as.numeric(Total$days),
                   PulACA = as.numeric(Total$Pulmonary.Adenocarcinoma),
                   HCC = as.numeric(Total$Hepatocellular.Carcinoma),
                   HSA = as.numeric(Total$Hemangiosarcoma),
                   HS = as.numeric(Total$Histiocytic.Sarcoma),
                   MammACA = as.numeric(Total$Mammary.Gland.Adenocarcinoma),
                   GCT = as.numeric(Total$Granulosa.Cell.Tumor),
                   Thyroid = as.numeric(Total$Thyroid.Tumor),
                   ThyroidAD = as.numeric(Total$Thyroid.Adenoma),
                   STS = as.numeric(Total$Soft.Tissue.Sarcomas),
                   AML = as.numeric(Total$Myeloid.Leukemia),
                   HardACA = as.numeric(Total$Harderian.Gland.Adenocarcinoma),
                   Harderian = as.numeric(Total$Harderian.Tumor),
                   HardAD = as.numeric(Total$Harderian.Gland.Adenoma),
                   LSA.BLL= as.numeric(Total$BLL),
                   LSA.Bmerge= as.numeric(Total$B.merge),
                   LSA.DLBCL= as.numeric(Total$DLBCL),
                   LSA.FBL= as.numeric(Total$FBL),
                   LSA.PreT = as.numeric(Total$PreT),
                   OSA = as.numeric(Total$Osteosarcoma),
                   PitAd = as.numeric(Total$Pituitary.Adenoma),
                   Amyloid = as.numeric(Total$Amyloidosis),
                   NN = as.numeric(Total$non.neoplastic),
                   ectoderm = as.numeric(Total$Ectoderm),
                   endoderm = as.numeric(Total$Endoderm),
                   mesoderm = as.numeric(Total$Mesoderm),
                   PSC = as.numeric(Total$Pulmonary.Sarcomatoid.Carcinoma),
                   Cat2 = as.numeric(Total$Cataract.2.0.Score.Event),
                   days2 = as.numeric(Total$Cataract.2.0.Score.Days))
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(row.names(pheno), "sex"))
HZE <- subset(pheno, group == "HZE")
Gamma <- subset(pheno, group == "Gamma")
Un <- subset(pheno, group == "Unirradiated")
All.irr <- subset(pheno, unirradiated == "0")



bootstrap <- HS.assoc.bootstrap(perms = 200, chr = 4, pheno = Gamma, pheno.col = "LSA.PreT",
                                probs, K, addcovar, markers, snp.file, outdir = "~/Desktop/files",
                                tx = "Gamma", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/",
                                peakMB = 82878119, window = 8000000)



Cat.Allirr.13 <- HS.cox.RMA.chrom(perms = 200, chr = 13, pheno = All.irr, pheno.col = "Cat2", days.col = "days2",
                           probs, K, addcovar, markers, snp.file, outdir = "~/Desktop/",
                           tx = "All Irradiated", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")


MammaACA.HZE <- HS.RMA.chrom(perms = 200, chr = 11, pheno = HZE, pheno.col = "MammACA", probs, K, addcovar,
                             markers, snp.file, outdir = "~/Desktop/",
                             tx = "HZE", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")
MammaACA.Gamma <- HS.RMA.chrom(perms = 200, chr = 11, pheno = Gamma, pheno.col = "MammACA", probs, K, addcovar,
                               markers, snp.file, outdir = "~/Desktop/",
                               tx = "Gamma", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")
MammaACA.Un <- HS.RMA.chrom(perms = 200, chr = 11, pheno = Un, pheno.col = "MammACA", probs, K, addcovar,
                               markers, snp.file, outdir = "~/Desktop/",
                               tx = "Unirradiated", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")
MammaACA.AI <- HS.RMA.chrom(perms = 200, chr = 11, pheno = All.irr, pheno.col = "MammACA", probs, K, addcovar,
                            markers, snp.file, outdir = "~/Desktop/",
                            tx = "All Irradiated", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")

### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
### PLOTTING BOOTSTRAP RESULTS (HISTOGRAM) ###
boot <- read.csv("~/Desktop/Whole Chr Resample/Raw Data-Cataract Chr 1.csv")
#Thy.boot2 = Thy.boot[which(Thy.boot$AML.LOD > 7),]

gamma = boot[which(boot$X == "gamma"),]
hze = boot[which(boot$X == "HZE"),]
un = boot[which(boot$X == "unirradiated"),]
#allirr = boot[which(Thy.boot$TX == "All.irradiated"),]

gamma = Thy.boot[which(Thy.boot$TX == "gamma" & Thy.boot$AML.LOD > 5),]
hze = Thy.boot[which(Thy.boot$TX == "HZE" & Thy.boot$AML.LOD > 5),]
un = Thy.boot[which(Thy.boot$TX == "Unirradiated" & Thy.boot$AML.LOD > 5),]
allirr = Thy.boot[which(Thy.boot$TX == "All.irradiated" & Thy.boot$AML.LOD > 5),]



#HISTOGRAM
layout(matrix(3:1, 3, 1))
hist(un$average, breaks=150, col="blue", main = "Unirradiated", xlab="Chromosome 1", prob = T, xlim = c(0, 125000000))
lines(density(un$average), col="black")
hist(gamma$average, breaks=150, col="green", main = "Gamma", xlab="Chromosome 1", prob = T, xlim = c(0, 125000000))
lines(density(gamma$average), col="black")
hist(hze$average, breaks=150, col="red", main = "HZE", xlab="Chromosome 1", prob = T, xlim = c(0, 125000000))
lines(density(hze$average), col="black", lwd = 1)

layout(matrix(4:1, 4, 1))
hist(un$average, breaks=150, col="blue", main = "Unirradiated", xlab="Chromosome 2", prob = T, xlim = c(0, 182113224))
lines(density(un$average), col="black")
hist(gamma$average, breaks=150, col="green", main = "Gamma", xlab="", prob = T, xlim = c(0, 182113224))
lines(density(gamma$average), col="black")
hist(hze$average, breaks=150, col="red", main = "HZE", xlab="", prob = T, xlim = c(0, 182113224))
lines(density(hze$average), col="black", lwd = 1)
hist(allirr$average, breaks=150, col="black", main = "All Irradiated", xlab="", prob = T, xlim = c(0, 182113224))
lines(density(All.irradiated$average), col="black", lwd = 1)

#HISTOGRAM WITHOUT DENSITY
layout(matrix(3:1, 3, 1))
hist(un$average, breaks=150, col="blue", main = "Unirradiated", xlab="Chromosome 2")
        #,xlim = c(110000000, 130000000))
hist(gamma$average, breaks=150, col="green", main = "Gamma", xlab="Chromosome 2")
        #,xlim = c(110000000, 130000000))
hist(hze$average, breaks=150, col="red", main = "HZE", xlab="Chromosome 2")
        #,xlim = c(110000000, 130000000))




#KERNAL DENSITY
layout(matrix(3:1, 3, 1))
d = density(hze$average)
plot(d, col="red", main = "HZE", xlab="Chromosome 2")
        #, xlim = c(0, 182113224))
d = density(gamma$average)
plot(d, col="green", main = "Gamma", xlab="Chromosome 2")
        #, xlim = c(0, 182113224))
d = density(un$average)
plot(d, col="blue", main = "Unirradiated", xlab="Chromosome 2")
        #, xlim = c(0, 182113224))

layout(matrix(4:1, 4, 1))
d = density(un$average)
plot(d, col="blue", main = "Unirradiated", xlab="Chromosome 2")
        #, xlim = c(0, 182113224))
d = density(gamma$average)
plot(d, col="green", main = "Gamma", xlab="Chromosome 2")
        #, xlim = c(0, 182113224))
d = density(hze$average)
plot(d, col="red", main = "HZE", xlab="Chromosome 2")
        #, xlim = c(0, 182113224))

d = density(allirr$average)
plot(d, col="blue", main = "Unirradiated", xlab="Chromosome 2", xlim = c(0, 182113224))

hist(hze$average, prob = T)
lines(density(x))


title(main = "Nonparametric Bootstrap Resampling with Replacement: Distribution of Peak LOD Scores")



boot <- read.csv("~/Desktop/Whole Chr Resample/Raw Data-MammACA 11.csv")

gamma = boot[which(boot$X == "gamma"),]
hze = boot[which(boot$X == "HZE"),]
un = boot[which(boot$X == "unirradiated"),]

merge.boot <- factor(boot$X, levels = c("gamma", "HZE", "unirradiated"),
                     labels = c("Gamma", "HZE", "Unirradiated"))
par(mfrow=c(1,1))
sm.density.compare(boot$average, boot$X, xlab = "Chromosome 11", lwd = 2.5)
title(main = "Resample Model Averaging: Mammary Adenocarcinoma")
colfill = c(2:(2+length(levels(merge.boot))))
legend("topleft", levels(merge.boot), fill = colfill)

layout(matrix(3:1, 3, 1))
hist(un$average, breaks=150, col="blue", main = "Unirradiated", xlab="Chromosome 11", prob = T, xlim = c(0, 125000000))
lines(density(un$average), col="black")
hist(gamma$average, breaks=150, col="green", main = "Gamma", xlab="", prob = T, xlim = c(0, 125000000))
lines(density(gamma$average), col="black")
hist(hze$average, breaks=150, col="red", main = "HZE", xlab="", prob = T, xlim = c(0, 125000000))
lines(density(hze$average), col="black", lwd = 1)




boot <- read.csv("~/Desktop/Whole Chr Resample/Raw Data-Cataract Chr 13.csv")
gamma = boot[which(boot$X == "gamma"),]
hze = boot[which(boot$X == "HZE"),]
un = boot[which(boot$X == "unirradiated"),]
allirr = boot[which(boot$X == "allirr"),]


merge.boot <- factor(boot$X, levels = c("allirr", "gamma", "HZE", "Unirradiated"),
                       labels = c("All Irradiated", "Gamma", "HZE", "Unirradiated"))
par(mfrow=c(1,1))
sm.density.compare(boot$average, boot$X, xlab = "Chromosome 13", lwd = 2.5)
title(main = "Resample Model Averaging: Cataractogenesis")
colfill = c(2:(2+length(levels(merge.boot))))
legend("topleft", levels(merge.boot), fill = colfill)



layout(matrix(4:1, 4, 1))
hist(un$average, breaks=150, col="blue", main = "Unirradiated", xlab="Chromosome 13", prob = T, xlim = c(0, 125000000))
lines(density(un$average), col="black")
hist(gamma$average, breaks=150, col="green", main = "Gamma", xlab="", prob = T, xlim = c(0, 125000000))
lines(density(gamma$average), col="black")
hist(hze$average, breaks=150, col="red", main = "HZE", xlab="", prob = T, xlim = c(0, 125000000))
lines(density(hze$average), col="black", lwd = 1)
hist(allirr$average, breaks=150, col="black", main = "All Irradiated", xlim = c(0, 125000000))
lines(density(allirr$average), col="black", lwd = 1)



