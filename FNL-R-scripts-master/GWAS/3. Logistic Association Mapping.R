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
load(file = "~/Desktop/R/QTL/WD/GRSD.Rdata")
load("/Users/elijah/Desktop/R/QTL/WD/hs.colors.Rdata")
setwd("~/Desktop/files")
outdir = "~/Desktop/files"

Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = Total$row.names, rownames = Total$row.names,
                   sex = as.numeric(Total$sex == "M"),
                   LG = as.numeric(Total$light.grey),
                   cohort = as.numeric(Total$Cohort),
                   group = as.character(Total$groups),
                   unirradiated = as.numeric(Total$Unirradiated),
                   days = as.numeric(Total$days),
                   AML = as.numeric(Total$Myeloid.Leukemia),
                   AML.ASXLdel = as.numeric(Total$Asxl1.Deletion),
                   AML.PU.1del = as.numeric(Total$Pu.1.Deletion))
                   
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
                   PulMet.transform = as.numeric(Total$PulMet.transform),
                   Metastatic.Tumors = as.numeric(Total$Metastatic.Tumors),
                   Pulmonary.Metastases = as.numeric(Total$Pulmonary.Metastases),
                   HCC.Metastatic.Density = as.numeric(Total$HCC.Metastatic.Density),
                   Tumors.that.could.met = as.numeric(Total$Tumors.that.could.met),
                   HCC...translocation = as.numeric(Total$HCC...translocation),
                   HCC.translocation = as.numeric(Total$HCC.translocation),
                   HCC.gel = as.numeric(Total$Gel.PCR),
                   cataract = as.numeric(Total$Cataract.2.0.Score.Event))
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(row.names(pheno), "sex"))

HZE <- subset(pheno, group == "HZE")
Gamma <- subset(pheno, group == "Gamma")
Un <- subset(pheno, group == "Unirradiated")
All.irr <- subset(pheno, unirradiated == "0")

Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = Total$row.names, rownames = Total$row.names,
                   sex = as.numeric(Total$sex == "M"),
                   group = Total$groups,
                   unirradiated = as.numeric(Total$Unirradiated),
                   days = as.numeric(Total$days),
                   HCC...translocation = as.numeric(Total$HCC...translocation),
                   HCC.translocation = as.numeric(Total$HCC.translocation),
                   HCC = as.numeric(Total$Hepatocellular.Carcinoma),
                   HCC.gel = as.numeric(Total$Gel.PCR))
#pheno = pheno[which(Total$Hepatocellular.Carcinoma=="1"),]
pheno = na.omit(pheno)
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(row.names(pheno), "sex"))
HZE <- subset(pheno, group == "HZE")
Gamma <- subset(pheno, group == "Gamma")
Un <- subset(pheno, group == "Unirradiated")
All.irr <- subset(pheno, unirradiated == "0")

HCC.met <- Total[which(Total$Hepatocellular.Carcinoma=="1" & Total$HCC.Metastatic.Density>0),]
HCC.met <- Total[which(Total$Hepatocellular.Carcinoma=="1"), ]
pheno = data.frame(row.names = HCC.met$row.names, sex = as.numeric(HCC.met$sex == "M"),
                   HCC.met = as.numeric(HCC.met$HCC.Met),
                   OSA = as.numeric(HCC.met$Osteosarcoma))

PulMET <- pheno[which(pheno$Tumors.that.could.met=="1"), ]

PSC <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/PSC.csv")
pheno = data.frame(row.names = PSC$row.names, rownames = PSC$row.names,
                   sex = as.numeric(PSC$Sex == "M"),
                   PSC = as.numeric(PSC$Sarcomatoid.Score))
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(row.names(pheno), "sex"))



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



