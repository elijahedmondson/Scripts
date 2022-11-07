library(DOQTL)

#Loading in data

load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))            

extract.raw.data(in.path = "/Volumes/External Terabyte/QTL/Founders", prefix = "",
                 out.path = "/Volumes/External Terabyte/QTL/extract.raw.data/Founders", 
                 array = "megamuga")

extract.raw.data(in.path = c("/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/1 - 96",
                             "/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/2 - 23 (last plate)",
                             "/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/2 - 481",
                             "/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/3 - 18 (last plate)",
                             "/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/4 - 600",
                             "/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/5 - 12 (last plate 600)",
                             "/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/6 - 648"), 
                 prefix = c("", "", "", "", "", "", ""),
                 out.path = "/Users/elijahedmondson/Desktop/R/QTL/extract.raw.data/GRSD", 
                 array = "megamuga")

#FIRST STEP

setwd("/Users/elijahedmondson/Desktop/R/QTL/WD")
getwd()
list.files("/Users/elijahedmondson/Desktop/R/QTL/WD")

# Read in founders.
fg = read.delim("/Users/elijahedmondson/Desktop/R/QTL/extract.raw.data/Founders/geno.txt")
fg[fg == TRUE] = "T"
fx = read.delim("/Users/elijahedmondson/Desktop/R/QTL/extract.raw.data/Founders/x.txt")
fy = read.delim("/Users/elijahedmondson/Desktop/R/QTL/extract.raw.data/Founders/y.txt")

# Load in data.
g = read.delim("/Users/elijahedmondson/Desktop/R/QTL/extract.raw.data/GRSD/geno.txt")
g[g == TRUE] = "T"
x = read.delim("/Users/elijahedmondson/Desktop/R/QTL/extract.raw.data/GRSD/x.txt")
y = read.delim("/Users/elijahedmondson/Desktop/R/QTL/extract.raw.data/GRSD/y.txt")

setwd("/Users/elijahedmondson/Desktop/R/QTL/WD/")
getwd()
save(x, y, g, fx, fy, fg, file = "GRSD1878data.Rdata")

load(file = "GRSD1878data.Rdata")

load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))

# Combine founders and samples.
x = rbind(fx, x)
y = rbind(fy, y)
g = rbind(fg, g)

# Remove outlier samples.
rmx = rowMeans(x, na.rm = T)
rmy = rowMeans(y, na.rm = T)
plot(rmx, rmy)
remove = which(rmx > 0.6)

View(remove)
list(remove)
[[1]]
97  237 1392  650 1614 1653 1697 1657 1671 1633 1707 1790 1636 
1639 1641 1686 1319 1429 1432 1414 1417 1419 1446 
627  783  903 1138 1245 1258 1259 1264 1282 1287 1289 1291 1293 
1299 1305 1306 1662 1794 1802 1825 1841 1849 1850 

x = x[-remove,]
y = y[-remove,]
g = g[-remove,]

sex = sex.predict(x, y, MM_snps, plot=T)
list(sex)
View(sex)


# All of the HS sample IDs are numbers.
fndr.idx = which(is.na(as.numeric(rownames(x))))
samp.idx = which(!is.na(as.numeric(rownames(x))))
fsex = sex[fndr.idx]
sex  = sex[samp.idx]

fx = x[fndr.idx,]
fy = y[fndr.idx,]
fg = g[fndr.idx,]
x = x[samp.idx,]
y = y[samp.idx,]
g = g[samp.idx,]

# A: A/J
# B: AKR/J
# C: BALB/cJ
# D: C3H/HeJ
# E: C57BL/6J
# F: CBA/J
# G: DBA/2J
# H: LP/J

#This needs to be specific for the data set and  in order
code = c("HH", "EE", "AA", "BB", "CC", "FF", "DD", "GG")  
names(code) = rownames(fx)
gen = rep(70, nrow(x))
names(gen) = rownames(x)

states = DOQTL:::create.genotype.states(LETTERS[1:8])

data = list(geno = as.matrix(g), sex = sex, gen = gen)
founders = list(geno = fg, sex = fsex, code = code, states = states)

# We only have male founders.
# For the allele call method, we're going to fake out the HMM by duplicating
# the males and calling them females.
founders$geno = as.matrix(rbind(founders$geno, founders$geno))
founders$sex = c(founders$sex, rep("F", length(founders$sex)))
names(founders$sex) = rownames(founders$geno)
founders$code = c(founders$code, founders$code)
names(founders$code) = rownames(founders$geno)

# 
attr(founders, "method") = "allele"
founders = add.missing.F1s(founders, MM_snps[,1:4]) 

getwd()
ls()

setwd("/Users/elijahedmondson/Desktop/R/QTL/WD")
getwd()
list.files("/Users/elijahedmondson/Desktop/R/QTL/WD")


calc.genoprob(data = data, chr = c(1:19, "X"), output.dir = "/Users/elijahedmondson/Desktop/R/QTL/HMM", 
              plot = T, array = "megamuga", sampletype = "HS", method = "allele", founders = founders)


plot.genoprobs(x = prsmth, snps = MM_snps, main = "1696")

recomb = summarize.genotype.transitions(path = "/Users/elijahedmondson/Desktop/R/QTL/HMM", snps = MM_snps)
head(recomb[[1]])
View(recomb)
mean(sapply(recomb, nrow))

#SECOND STEP
library(DOQTL)

setwd("/Users/elijahedmondson/Desktop/R/QTL/WD")
getwd()
list.files("/Users/elijahedmondson/Desktop/R/QTL/WD")

load(file = "/Users/elijahedmondson/Desktop/R/QTL/HMM/founder.probs.Rdata")
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))


#EFE add "bychr = TRUE"
K = kinship.probs(model.probs, bychr = TRUE, snps = MM_snps)
#can you plot?

pheno = data.frame(row.names = Sarcomatoid$row.names, sex = as.numeric(Sarcomatoid$Sex == "M"),  
                   sarcomatoid = as.numeric(Sarcomatoid$sarcomatoid))

pheno = data.frame(row.names = GRSD.phenotype$row.names, sex = as.numeric(GRSD.phenotype$Sex == "M"), 
                   Age = as.numeric(GRSD.phenotype$days), 
                   Weight = GRSD.phenotype$weight, 
                   BCS = GRSD.phenotype$BCS, 
                   Albino = as.numeric(GRSD.phenotype$albino), 
                   Amyloidosis = as.numeric(GRSD.phenotype$amyloidosis),
                   Lymphoma = as.numeric(GRSD.phenotype$lymphoma), 
                   Harderian = as.numeric(GRSD.phenotype$harderian.gland),
                   Histiocytic.sarcoma = as.numeric(GRSD.phenotype$histiocytic.sarcoma),
                   AML = as.numeric(GRSD.phenotype$myeloid.leukemia),
                   LSA.AML = as.numeric(GRSD.phenotype$lymphoma.or.leukemia),
                   Lymphoid.leukemia = as.numeric(GRSD.phenotype$lymphoid.leukemia),
                   Pituitary.adenoma = as.numeric(GRSD.phenotype$pituitary.adenoma),
                   HCC = as.numeric(GRSD.phenotype$HCC),
                   Pulmonary.adenocarcinoma = as.numeric(GRSD.phenotype$pulmonary.adenocarcinoma),
                   Mammary.adenocarcinoma = as.numeric(GRSD.phenotype$mammary.adenocarcinoma),
                   GCT = as.numeric(GRSD.phenotype$GCT),
                   OSA = as.numeric(GRSD.phenotype$OSA),
                   Sarcoma = as.numeric(GRSD.phenotype$sarcoma),
                   HSA = as.numeric(GRSD.phenotype$HSA),
                   Metastasis = as.numeric(GRSD.phenotype$metastasis),
                   Solid.tumors = as.numeric(GRSD.phenotype$Solid.tumors))
                   

covar = data.frame(sex = as.numeric(Sarcomatoid$Sex == "M"))
rownames(covar) = rownames(pheno)

qtl = scanone(pheno = pheno, pheno.col = "sarcomatoid", probs = model.probs, K = K, 
              addcovar = covar, snps = MM_snps)

save(qtl, file = "")

plot(qtl, main = "EMT")

perms = scanone.perm(pheno = pheno, pheno.col = "sarcomatoid", probs = model.probs, addcovar = covar, 
                     snps = MM_snps, path = "/Users/elijahedmondson/Desktop/R/QTL/perms", 
                     nperm = 1000)

save(perms, file = "Sarcomatoid.Rdata")
#scanone.perm(n=10) takes about 1 hour

thr1 = quantile(perms, probs = 0.90)
thr2 = quantile(perms, probs = 0.95)
thr3 = quantile(perms, probs = 0.99)

plot(AM.qtl, sig.thr = c(thr1, thr2, thr3), main = "Epithelial to Mesenchymal Transition")

setwd("/Users/elijahedmondson/Desktop/R/QTL/WD")
getwd()
save(perms, file = "SarcomatoidPerms.Rdata")
save(qtl, file = "SarcomatoidQTL.Rdata")
list.files("/Users/elijahedmondson/Desktop/R/QTL/WD")

interval = bayesint(AM.qtl, chr = 14)
interval
mgi = get.mgi.features(chr = interval[1,2], start = interval[1,3], end = interval[3,3], 
                       type = "gene", source = "MGI")
nrow(mgi)
head(mgi)


ma = assoc.map(pheno = pheno, pheno.col = "albino", probs = model.probs, K = K, addcovar = covar, 
               snps = MM_snps, chr = 7, start = 87, end = 88)
coefplot(qtl, chr = 1)
tmp = assoc.plot(ma, thr = 7)
unique(tmp$sdps)


load(file = "/Users/elijahedmondson/Desktop/R/QTL/WD/GRSD.K.model.probs.RData")
load(file = "")
load(file = "")
load(file = "")
library(DOQTL)




