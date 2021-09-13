# LOAD PACKAGES #
library(BSgenome.Mmusculus.UCSC.mm10)
library(doParallel)
library(VariantAnnotation)
library(GenomicRanges)
library(regress)
library(MASS)
library(DOQTL)
library(lmtest)
library(HZE)
library(dplyr)
options(stringsAsFactors = F)
C:\Users\edmondsonef\Downloads

load(file = "C:/Users/edmondsonef/Downloads/HZE.new.Rdata")
load("/Users/elijah/Desktop/R/QTL/WD/hs.colors.Rdata")
setwd("~/Desktop/files")
outdir = "~/Desktop/files"



GRSD.assoc(pheno = Un, pheno.col = "cataract", probs, K, addcovar = addcovar,
           markers, snp.file = "snp.file", outdir = "~/Desktop/files", tx = "Unirradiated",
           sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")

GRSD.poisson(pheno = Un, pheno.col = "HCC...translocation", probs, K, addcovar = addcovar,
             markers, snp.file = "snp.file", outdir = "~/Desktop/files", tx = "Unirradiated",
             sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")


perms <- GRSDassoc.perms(perms = 2, chr = 19, pheno = HZE, Xchr = F, addcovar = addcovar,
                         pheno.col = "HCC", probs = probs, K = K, markers = markers,
                         snp.file = snp.file, outdir = "~/Desktop/files", tx = "Test",
                         sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")

bootstrap <- HS.assoc.bootstrap(perms = 200, chr = 3, pheno = pheno, pheno.col = "HCC...translocation",
                                probs, K, addcovar, markers, snp.file, outdir = "~/Desktop/files",
                                tx = "ALL", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/",
                                peakMB = 55139932)
