
library(BSgenome.Mmusculus.UCSC.mm10)
library(doParallel)
library(foreach)
library(Rsamtools)
library(VariantAnnotation)
library(DOQTL)
library(GenomicRanges)
library(regress)
library(MASS)
library(lmtest)
library(HZE)
options(stringsAsFactors = F)
load(file = "/home/ubuntu/CloudLOG.Rdata")
setwd("~/Dropbox/Rstudio.cloud/WD")
outdir = "~/Dropbox/Rstudio.cloud/WD"


perms <- GRSDassoc.perms(perms = 2, chr = 1:2, pheno = HZE, Xchr = F, addcovar = HZE.add,
                         pheno.col = "HCC", probs = probs, K = K, markers = markers,
                         snp.file = snp.file, outdir = "~/Dropbox/Rstudio.cloud/WD", tx = "Test",
                         sanger.dir = "/home/ubuntu/HS.sanger.files/")


source("https://bioconductor.org/biocLite.R")
biocLite("doParallel", "foreach", "Rsamtools", "VariantAnnotation", "GenomicRanges")
