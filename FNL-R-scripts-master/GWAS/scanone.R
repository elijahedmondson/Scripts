#SECOND STEP
library(DOQTL)

setwd("/Users/elijahedmondson/Desktop/R/QTL/WD")
getwd()
list.files("/Users/elijahedmondson/Desktop/R/QTL/WD")

load(file = "/Users/elijahedmondson/Desktop/R/QTL/HMM/founder.probs.Rdata")
load(file = "/Users/elijahedmondson/Desktop/R/QTL/WD/kinship.probs1837.RData")
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))


pheno = data.frame(row.names = Sheet.1.Table.1$row.names, sex = as.numeric(Sheet.1.Table.1$Sex == "M"), 
                   osteosarcoma = as.numeric(Sheet.1.Table.1$osteosarcoma))

covar = data.frame(sex = as.numeric(Sheet.1.Table.1$Sex == "M"))
rownames(covar) = rownames(pheno)

pheno = data.frame(row.names = GRSD.phenotype$row.names, sex = as.numeric(GRSD.phenotype$Sex == "M"), 
                   age = as.numeric(GRSD.phenotype$age.in.days), 
                   weight = as.numeric(GRSD.phenotype$weight), 
                   black = GRSD.phenotype$black, 
                   albino = as.numeric(GRSD.phenotype$albino),
                   lymphoma = as.numeric(GRSD.phenotype$Thymic.Lymphoma),
                   HCC = as.numeric(GRSD.phenotype$HCC))

covar = data.frame(sex = as.numeric(GRSD.phenotype$Sex == "M"))
rownames(covar) = rownames(pheno)

qtl = scanone(pheno = pheno, pheno.col = "osteosarcoma", probs = model.probs, K = K, 
              addcovar = covar, snps = MM_snps)

save(qtl, file = "OSAqtl.Rdata")

plot(qtl, main = "Osteosarcoma")

perms = scanone.perm(pheno = pheno, pheno.col = "osteosarcoma", probs = model.probs, addcovar = covar, 
                     snps = MM_snps, path = "/Users/elijahedmondson/Desktop/R/QTL/perms", nperm = 500)

save(perms, file = "OSAPerms.Rdata")

thr = quantile(perms, probs = 0.9)
thr2 = quantile(perms, probs = 0.95)
thr3 = quantile(perms, probs = 0.99)

plot(qtl, sig.thr = c(thr2, thr, thr3), main = "OSA")

coefplot(qtl, chr = 14)
interval = bayesint(qtl, chr = 14)
interval
ma = assoc.map(pheno = pheno, pheno.col = "black", probs = model.probs, K = K, addcovar = covar, 
               snps = MM_snps, chr = interval[1,2], start = interval[1,3], end = interval[3,3])

tmp = assoc.plot(ma, thr = 11)
unique(tmp$sdps)
mgi = get.mgi.features(chr = interval[1,2], start = interval[1,3], end = interval[3,3], type = "gene", source = "MGI")
nrow(mgi)
head(mgi)

save(mgi, file = "Age.mgi.Rdata")
