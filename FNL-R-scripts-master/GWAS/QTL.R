library(DOQTL)

load(file = "~/Desktop/R/Build/K.Rdata")
load(file = "~/Desktop/R/Build/model.probs.Rdata")
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))

save(MM_snps, K, model.probs, sdp.file, file = "GRSD_master.Rdata")
load(file = "/Users/elijahedmondson/Desktop/R/QTL/WD/GRSD_master.Rdata")



Sarcomatoid = read.csv("/Users/elijahedmondson/Desktop/R/GRSD.phenotype/CSV/PSC.csv")

pheno = data.frame(row.names = Sarcomatoid$row.names, sex = as.numeric(Sarcomatoid$Sex == "M"),  
                   sarcomatoid = as.numeric(Sarcomatoid$Sarcomatoid.Score))

covar = data.frame(sex = as.numeric(Sarcomatoid$Sex == "M"))
rownames(covar) = rownames(pheno)

GRSD.phenotype = read.csv("/Users/elijahedmondson/Desktop/R/GRSD.phenotype/CSV/GRSD.phenotype.csv")

pheno = data.frame(row.names = GRSD.phenotype$row.names, sex = as.numeric(GRSD.phenotype$sex == "M"), 
                   Age = as.numeric(GRSD.phenotype$days), 
                   Weight = GRSD.phenotype$weight, 
                   BCS = GRSD.phenotype$BCS, 
                   Albino = as.numeric(GRSD.phenotype$albino),
                   Black = as.numeric(GRSD.phenotype$black),
                   AML = as.numeric(GRSD.phenotype$Myeloid.Leukemia),
                   Lymphoma = as.numeric(GRSD.phenotype$Lymphoma))

covar = data.frame(sex = as.numeric(GRSD.phenotype$sex == "M"))
rownames(covar) = rownames(pheno)

qtl = scanone(pheno = pheno, pheno.col = "Lymphoma", probs = model.probs, K = K, 
              addcovar = covar, snps = MM_snps)

plot(qtl, main = "")

perms = scanone.perm(pheno = pheno, pheno.col = "Amyloidosis", probs = model.probs, addcovar = covar, 
                     snps = MM_snps, path = "/Users/elijahedmondson/Desktop/R/QTL/perms", 
                     nperm = 1000)

thr1 = quantile(perms, probs = 0.90)
thr2 = quantile(perms, probs = 0.95)
thr3 = quantile(perms, probs = 0.99)

plot(qtl, sig.thr = c(thr1, thr2, thr3), main = "")

setwd("/Users/elijahedmondson/Desktop/R/QTL/WD")
getwd()
save(perms, file = "OSAPerms.Rdata")
save(qtl, file = "LSA.AMLQTL.Rdata")
list.files("/Users/elijahedmondson/Desktop/R/QTL/WD")

interval = bayesint(qtl, chr = 14)
interval
mgi = get.mgi.features(chr = interval[1,2], start = interval[1,3], end = interval[3,3], type = "gene", source = "MGI")
nrow(mgi)
head(mgi)


ma = assoc.map(pheno = pheno, pheno.col = "sarcomatoid", probs = model.probs, K = K, addcovar = covar, 
               snps = MM_snps, chr = interval[1,2], start = interval[1,3], end = interval[3,3])
coefplot(qtl, chr = 1)
tmp = assoc.plot(ma, thr = 1)
unique(tmp$sdps)


load(file = "/Users/elijahedmondson/Desktop/R/QTL/WD/GRSD.K.model.probs.RData")
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
load(file = "/Users/elijahedmondson/Desktop/R/QTL/HMM/chr19.emission.probs.Rdata")
load(file = "")
library(DOQTL)

# ASSOCIATION FROM D GATTI 2/4/15 #

library(DOQTL)
setwd("/Users/elijahedmondson/Desktop/R/QTL/WD")

sdp.file = "/Users/elijahedmondson/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz"

pheno[,2] = as.numeric(pheno[,2]) - 1

qtlscan = scanone.assoc(pheno = pheno, pheno.col = "", probs = model.probs, K = K, 
                    addcovar = covar, markers = MM_snps, sdp.file = sdp.file, ncl = 2)

save(qtlscan, file = "Amyloidosis_AM.Rdata")

png("QTL.png", width = 2000, height = 1600, res = 128)
DOQTL:::plot.scanone.assoc(qtlscan, bin.size = 100)
dev.off()

png("sarcomatoid_QTL_chr1.png", width = 2000, height = 1600, res = 128)
DOQTL:::plot.scanone.assoc(qtlscan, chr = 6, bin.size = 10)
dev.off()


