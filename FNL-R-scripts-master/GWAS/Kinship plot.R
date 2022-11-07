

library(BSgenome.Mmusculus.UCSC.mm10)
library(doParallel)
library(foreach)
library(Rsamtools)
library(DOQTL)
library(VariantAnnotation)
library(GenomicRanges)
library(regress)
library(MASS)
library(lmtest)
library(HZE)
library(reshape2)
library(ggplot2)


load(file = "~/Desktop/R/QTL/WD/HS\ HMM\ Rdata/K.Rdata")
load(file = "~/Desktop/R/QTL/WD/HS\ HMM\ Rdata/model.probs.Rdata")

csv.m <- melt(K6, id.vars="V1")
csv.m$V1 <- factor(csv.m$V1, levels=unique(as.character(csv.m$V1)) )

qplot(x=Var2, y=Var1, data=csv.m, fill=value, geom="tile") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
        scale_fill_gradient2(low = "blue", mid = "green", high = "red")

K6 <- K[6]


