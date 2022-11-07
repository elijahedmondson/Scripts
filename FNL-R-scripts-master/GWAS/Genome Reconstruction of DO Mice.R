library(DOQTL)

'y' <- read.delim("/Users/elijahedmondson/Desktop/R/QTL/extract.raw.data/GRSD/y.txt", header = T)
'x' <- read.delim("/Users/elijahedmondson/Desktop/R/QTL/extract.raw.data/GRSD/x.txt", header = T)
geno <- read.delim("/Users/elijahedmondson/Desktop/R/QTL/extract.raw.data/GRSD/geno.txt", header = T)
pheno <- read.delim("/Users/elijahedmondson/Desktop/R/GRSD.phenotype/CSV/GRSD.csv", header = T, sep = ",")

load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))

sex = sex.predict(x = x, y = y, snps = MM_snps, plot = T)

gen = paste("DO", gsub("70", "", pheno$Gen), sep = "")
names(gen) = rownames(pheno)
gen = gen[names(gen) %in% names(sex)]
gen = gen[match(names(sex), names(gen))]
stopifnot(all(rownames(x) == names(sex)))
stopifnot(all(rownames(x) == names(gen)))
data = list(x = x, y = y, sex = sex, gen = gen)

#ALL BELOW IS FOR "FOUNDERS"
'yf' <- read.delim("~/Desktop/R/CSU Geneseek/Founders.raw.data/y.txt", header = T)
'xf' <- read.delim("~/Desktop/R/CSU Geneseek/Founders.raw.data/x.txt", header = T)
genof <-  read.delim("~/Desktop/R/CSU Geneseek/Founders.raw.data/geno.txt", header = T)
  
sexf = sex.predict(x = xf, y = yf, snps = MM_snps, plot = T)  

founders = list(x = xf, y = yf, sex = sexf)
#ALL ABOVE IS FOR FOUNDERS


rm(x, y, sex, gen)
gc()
calc.genoprob(data = data, chr = c(1:19, "X"), 
              output.dir = "~/Desktop/R/CSU Geneseek/NEW/", 
              plot = T, array = "megamuga", sampletype = "CC", 
              method = "intensity", founders = founders, transprobs, snps)

founders = list(x = xf, y = yf)
