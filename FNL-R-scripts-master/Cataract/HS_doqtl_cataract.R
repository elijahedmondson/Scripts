library(BSgenome.Mmusculus.UCSC.mm10)
library(doParallel)
library(foreach)
library(Rsamtools)
library(VariantAnnotation)
library(DOQTL)
library(GenomicRanges)
library(survival)
library(regress)
library(HZE)

load(file = "C:/Users/edmondsonef/Desktop/QTL/HZEproject.Rdata")
load(file = "C:/Users/edmondsonef/Desktop/QTL/HS.colors.Rdata")
#load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
sdp.file = "C:/Users/edmondsonef/Desktop/QTL/HS_Sanger_SDPs.txt.bgz"

outdir = "C:/Users/edmondsonef/Desktop/R-plots/"

total <- read_excel("C:/Users/edmondsonef/Desktop/Cataract/CATARACT_final.xlsx", sheet ="simple")
pheno <- data.frame(row.names = total$row.names,
                    sex = as.numeric(total$sex == "M"),
                    group = as.character(total$groups),
                    albino = as.numeric(total$`coat color` == "albino"),
                    cat_average = as.numeric(total$`cat_average`),
                    cat_sum = as.numeric(total$`cat_sum`),
                    cat_difference = as.numeric(total$`cat_difference`))
