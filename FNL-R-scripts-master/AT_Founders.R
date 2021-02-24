
#Heterogeneous stock (HS)
#The cross type for heterogenous stock (HS) derived from 8 inbred strains is "hs". HS is similar to AIL, but starting with 8 inbred strains 
#rather than 2.
#
#For HS and other multi-parent populations, we need SNP genotypes for both the individuals and the founders. They'll be coded as 1, 2, and 3, 
#with 1 and 3 being the two homozygotes and 2 the heterozygote. In the genotype and founder genotype files, the genotypes need to have a 
#consistent coding scheme (like A, H, and B) for all markers. While the data might arrive with nucleotide codes, they need to be recoded before 
#they can be read into R. (The R/qtl2convert package includes some tools to assist with this.)
#
#The yaml control file will includes lines similar to those for the F2 intercross:
#  
# genotypes:
#A: 1
#H: 2
#B: 3
#The json equivalent is the following:
#  
#  "genotypes": {"A":1, "H": 2, "B":3}
#The genotype codes need not be single-character codes; for example, you could use AA, AB, and BB.
#
#A single cross info covariate is needed: the generation number for each individual. This could be included with the other covariates, or it can be in a separate file.
#
#If the X chromosome is included, a sex covariate is also needed.









founder <- read.delim("C:/Users/edmondsonef/Desktop/AT_HS/AT_HS.NPT_Snps/AT_SNPs/CSU_batch1/batch1_FinalReport/batch1_FinalReport_noheader.txt")

batch2 <- read.delim("C:/Users/edmondsonef/Desktop/AT_HS/AT_HS.NPT_Snps/AT_SNPs/CSU_batch2/batch2_FinalReport/batch2_FinalReport_noheader.txt")

# subset founder genotypes

founders <- subset(batch2, batch2$Sample.ID %in% c('A/J F', 'A/J M', 'C57 F', 'C57 M', 'CBA M', 'CBA F',
                                               'C3H F', 'C3H M', 'BALB F', 'BALB M', 'DBA M', 'DBA F',
                                               'AKR F', 'AKR M', 'LP M', 'LP F'))
Founders <- rbind(founder, founders)

write.table(Founders, file = "C:/Users/edmondsonef/Desktop/AT_HS/AT_HS.NPT_Snps/AT_SNPs/CSU_batch1/batch1_FinalReport/Founders.8xMF.LPxM.17.txt", 
            quote = FALSE, row.names = FALSE, sep = "\t")

######EFE#######



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
