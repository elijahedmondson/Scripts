library(devtools)
load_all("d:/182_DO_QTL_Mapping/DOQTL")

setwd("C:/Users/dgatti/Documents/HS/")

# Read in founders.
#fg = read.delim("Founders/geno.txt")
#fx = read.delim("Founders/x.txt")
#fy = read.delim("Founders/y.txt")

# Load in data.
#g = read.delim("Samples/geno.txt")
#x = read.delim("Samples/x.txt")
#y = read.delim("Samples/y.txt")

#save(x, y, g, fx, fy, fg, file = "HSdata.Rata")

load(file = "HSdata.Rata")

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
x = x[-remove,]
y = y[-remove,]
g = g[-remove,]

sex = sex.predict(x, y, MM_snps, T)

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

code = c("HH", "EE", "AA", "BB", "CC", "FF", "DD", "GG", "AA", "BB", "CC")
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

save(founders, data, file = "HS_allele.Rdata")

load(file = "HS_allele.Rdata")

setwd("HMM")
calc.genoprob(data = data, chr = 1:19, output.dir = ".", plot = T, array = "megamuga",
              sampletype = "HS", method = "allele", founders = founders)

