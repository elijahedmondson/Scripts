# Lifespan, females only, COXPH with unique SDPs.
library(DOQTL)
library(doParallel)
library(foreach)
library(Rsamtools)
library(VariantAnnotation)
library(GenomicRanges)
library(survival)
library(regress)
options(stringsAsFactors = F)

### DMG
#setwd("/Users/elijahedmondson/Desktop/R/QTL/WD")
setwd("/hpcdata/dgatti/HS/fromElijah/")

# Pass in the number of clusters (nodes).
#args = commandArgs(trailingOnly = TRUE)
#ncl = as.numeric(args[[1]])

### DMG
#ncl = 2
ncl = 4

####################
# THINGS TO CHANGE

# Set the output file directory.
###DMG
#outdir = "/Users/elijahedmondson/Desktop/R/QTL/WD/hq_snps"
outdir = "/hpcdata/dgatti/HS/fromElijah/QTL"

# Load in your data. This file contains pheno, probs, markers and K.
###DMG
#HZE <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/HZE.csv")
#load(file = "/Users/elijahedmondson/Desktop/R/QTL/WD/model.probs.Rdata")
#probs <- model.probs
#rm(model.probs)

#load(file = "/Users/elijahedmondson/Desktop/R/QTL/WD/K.Rdata")

#pheno = data.frame(row.names = HZE$row.names, sex = as.numeric(HZE$sex == "M"),  
#                   days = as.numeric(HZE$Cataract.2.0.Score),
#                   cataract = as.numeric(HZE$Cataract.2.0.Event),
#                   LSA = as.numeric(HZE$Lymphoma))

load("HZE.Rdata")

###DMG: you have sex coded as 0/1, so this needs to change.
### I set the covariate below.
#covar = data.frame(sex = as.numeric(pheno$sex == "M"))
#addcovar = cbind(sex = as.numeric(factor(pheno$sex)) - 1)
#rownames(addcovar) = rownames(pheno)
#rownames(covar) = rownames(pheno)
#rm(HZE)

load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
markers <- MM_snps
rm(MM_snps)
markers = markers[markers[,1] %in% dimnames(probs)[[3]],]
###DMG
### I had the marker positions in Mb and the Sanger SNP positions
### in bp below.
markers[,3] = markers[,3] * 1e6

stopifnot(nrow(markers) == dim(probs)[3])
stopifnot(markers[,1] == dimnames(probs)[[3]])

# Set the Sanger SNP file location.
###DMG
### My fault. We don't use the SDP file here. Sorry, I forgot which
### version of the code you had.
#sdp.file = "/Users/elijahedmondson/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz"
snp.file = "ftp://ftp.jax.org/SNPtools/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"

# Set the cross type (DO or HS).
cross = "HS"

# Set a file prefix for the output files.
file.prefix = "Cataract_Latency_coxph"

# Set the plot title.
plot.title = "Cataract Latency, CoxPH, HQ SNPs"

###DMG
### I moved this out of the workfxn.
### CHANGE THIS FOR EACH SURVIVAL PHENOTYPE.
# Create the survival function.
surv = Surv(pheno$days, pheno$cataract)

####################

# Set up covariates
addcovar = matrix(pheno$sex, ncol = 1, dimnames = 
          list(rownames(pheno), "sex"))

###DMG
## Remove samples that are not found in both pheno and probs.
#pheno = pheno[rownames(pheno) %in% dimnames(probs)[[1]],,drop = FALSE]
#probs = probs[dimnames(probs)[[1]] %in% rownames(pheno),,]
#probs = probs[match(rownames(pheno), dimnames(probs)[[1]]),,]
#probs = probs[,,dimnames(probs)[[3]] %in% markers[,1]]
samples = intersect(rownames(pheno), rownames(probs))
samples = intersect(samples, rownames(addcovar))
samples = intersect(samples, rownames(K[[1]]))

stopifnot(length(samples) > 0)
print(paste("Found", length(samples), "samples in common."))

# Synch up the samples in all of the data.
pheno = pheno[samples,,drop = FALSE]
addcovar = addcovar[samples,,drop = FALSE]
probs = probs[samples,,,drop = FALSE]
for(i in 1:length(K)) {
  K[[i]] = K[[i]][samples, samples]
} # for(i)

# Split up the data by chromosome.
chrs = c(1:19, "X")
data = vector("list", length(chrs))
names(data) = chrs
for(i in 1:length(chrs)) {

  rng = which(markers[,2] == chrs[i])
  data[[i]] = list(probs = probs[,,rng], K = K[[i]],
              markers = markers[rng,])

} # for(i)
rm(probs, K, markers)

# Make a plot of the survival data.
fit = survfit(surv ~ addcovar)
plot(fit, col = 1:2, las = 1, main = plot.title)
legend("bottomleft", col = 1:2, lty = 1, legend = c("female", "male"))
mod = coxph(surv ~ addcovar)
text(x = 25, y = 0.15, labels = paste("p =", format(anova(mod)[2,4],
     digits = 2)), adj = 0)

setwd(outdir)

# Make a function for each worker to execute.
### DMG
### I changed the name of the argument to 'obj' so that 
### you don't confuse it with the 'data' list.
workfxn = function(obj) {

  chr = obj$markers[1,2]

  setwd(outdir)

  # Get the Sanger SNPs.
###DMG
### You're working with the HS, so just get HS colors.
#  strains = sub("/", "_", do.colors[,2])
#  if(cross == "HS") {
    strains = sub("/", "_", hs.colors[,2])
#  } # if(cross = "HS")

  # Read the Sanger VCF file.
  hdr = scanVcfHeader(snp.file)
  gr = GRanges(seqnames = chr, range = IRanges(start = 0, 
       end = 200e6))
  param = ScanVcfParam(geno = c("GT", "FI"), fixed = "ALT",
          samples = strains[strains != "C57BL_6J"], which = gr)
  sanger = readVcf(file = snp.file, genome = "mm10", param = param)

  # Keep high quality SNPs (quality == 1)
  sanger = sanger[rowSums(geno(sanger)$FI, na.rm = TRUE) == 7]

  # Keep polymorphic SNPs.
  keep = which(rowSums(geno(sanger)$GT == "0/0", na.rm = TRUE) < 7)
  sanger = sanger[keep]
  rm(keep)

  # We have to do some work to extract the alternate allele.
  alt = CharacterList(fixed(sanger)$ALT)
  alt = unstrsplit(alt, sep = ",")

  # Extract the SNP positions and genotypes.
###DMG
### Changed 'rowData()' to 'rowRanges()' because rowData was deprecated.
  sanger.hdr = data.frame(ID = names(rowRanges(sanger)), CHR = as.character(seqnames(sanger)),
               POS = start(sanger), REF = as.character(fixed(sanger)$REF),
               ALT = alt, stringsAsFactors = FALSE)
  rm(alt)

###DMG
### Again, you have HS mice. Just use the HS data. You can delete the DO lines.
  # Add C57BL/6J to the Sanger SNPs.
#  if(cross == "DO") {
#    sanger = cbind("A_J" = geno(sanger)$GT[,1,drop = FALSE], 
#             "C57BL_6J" = "0/0",
#             geno(sanger)$GT[,2:7,drop = FALSE])
#  } else if(cross == "HS") {
    sanger = cbind(geno(sanger)$GT[,1:4,drop = FALSE],
             "C57BL_6J" = "0/0",
             geno(sanger)$GT[,5:7,drop = FALSE])
#  } # else

  # Convert allele calls to numeric values.
  sanger = (sanger != "0/0") * 1

  # Make the MAF between 1/8 and 4/8.
  flip = which(rowSums(sanger) > 4)
  sanger[flip,] = 1 - sanger[flip,,drop = FALSE]
  rm(flip)

###DMG
### I'm moving this outside of the function.
  # Create the survival object.
#  surv = Surv(pheno$days, pheno$cataract)

  # Null model.
###DMG
### Put the null logistic regression or linear model here.
  null.mod = coxph(surv ~ addcovar)
  null.ll = logLik(null.mod)
  pv = rep(0, nrow(sanger))

  # Get the unique SDPs between each pair of markers and
  # calculate the COXPH LOD.

  # CoxPH function.
  coxph.fxn = function(snp.rng, local.probs) {

    # Get the SDPs.
    sdp.nums = sanger[snp.rng,] %*% 2^(7:0)
    sdps2keep = which(!duplicated(sdp.nums))
    cur.sdps = sanger[snp.rng,,drop = FALSE][sdps2keep,,drop = FALSE]
    unique.sdp.nums = sdp.nums[sdps2keep]
    m = match(sdp.nums, unique.sdp.nums)

    # Multiply the SDPs by the haplotype probabilities.
    cur.alleles = tcrossprod(cur.sdps, local.probs)
    cur.ll = rep(null.ll, nrow(cur.sdps))

    # Check for low allele frequencies and remove SDPs with too
    # few samples carrying one allele.
    sdps.to.use = which(rowSums(cur.alleles) > 0.5)

    # Run the Cox PH model at each unique SDP.
    for(j in sdps.to.use) {

###DMG
### Put the logistic regression or linear model here.

      mod = coxph(surv ~ addcovar + cur.alleles[j,])
      cur.ll[j] = logLik(mod)

    } # for(j)

    # This is the LRS.
    cur.ll = cur.ll - null.ll

    # Return the results.
    cur.ll[m]

   } # coxph.fxn()

   # SNPs before the first marker.
   snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
   if(length(snp.rng) > 0) {

     pv[snp.rng] = coxph.fxn(snp.rng, obj$probs[,,1])

   } # if(length(snp.rng) > 0)

   # SNPs between Markers.
   for(i in 1:(nrow(obj$markers)-1)) {

     snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                     sanger.hdr$POS <= obj$markers[i+1,3])

     if(length(snp.rng) > 0) {

       # Take the mean of the haplotype probs at the surrounding markers.
       pv[snp.rng] = coxph.fxn(snp.rng, (obj$probs[,,i] + 
                     obj$probs[,,i+1]) * 0.5)

     } # if(length(snp.rng) > 0)

   } # for(i)

  # SNPs after the last marker.
  snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
  if(length(snp.rng) > 0) {

    pv[snp.rng] = coxph.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])

  } # if(length(snp.rng) > 0)

  # Convert LRS to p-values using the chi-squared distribution.
  pv = pchisq(2 * pv, df = 7, lower.tail = FALSE)
  pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)

  save(pv, file = paste0(file.prefix, "_chr", chr, ".Rdata"))

  png(paste0(file.prefix, "_chr", chr,".png"), width = 2000, 
      height = 1600, res = 200)
  plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
  mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
  dev.off()

  # Return the positions and p-values.
  return(pv)

} # workfxn()


### DMG
# Special function to map the X chromosome correctly.
# We map using sex as an interactive covariate.
workfxn.xchr = function(obj) {

  chr = obj$markers[1,2]

  setwd(outdir)

  # Get the Sanger SNPs.
###DMG
### You're working with the HS, so just get HS colors.
#  strains = sub("/", "_", do.colors[,2])
#  if(cross == "HS") {
    strains = sub("/", "_", hs.colors[,2])
#  } # if(cross = "HS")

  # Read the Sanger VCF file.
  hdr = scanVcfHeader(snp.file)
  gr = GRanges(seqnames = chr, range = IRanges(start = 0, 
       end = 200e6))
  param = ScanVcfParam(geno = c("GT", "FI"), fixed = "ALT",
          samples = strains[strains != "C57BL_6J"], which = gr)
  sanger = readVcf(file = snp.file, genome = "mm10", param = param)

  # Keep high quality SNPs (quality == 1)
  sanger = sanger[rowSums(geno(sanger)$FI, na.rm = TRUE) == 7]

  # Keep polymorphic SNPs.
  keep = which(rowSums(geno(sanger)$GT == "0/0", na.rm = TRUE) < 7)
  sanger = sanger[keep]
  rm(keep)

  # We have to do some work to extract the alternate allele.
  alt = CharacterList(fixed(sanger)$ALT)
  alt = unstrsplit(alt, sep = ",")

  # Extract the SNP positions and genotypes.
###DMG
### Changed 'rowData()' to 'rowRanges()' because rowData was deprecated.
  sanger.hdr = data.frame(ID = names(rowRanges(sanger)), CHR = as.character(seqnames(sanger)),
               POS = start(sanger), REF = as.character(fixed(sanger)$REF),
               ALT = alt, stringsAsFactors = FALSE)
  rm(alt)

###DMG
### Again, you have HS mice. Just use the HS data. You can delete the DO lines.
  # Add C57BL/6J to the Sanger SNPs.
#  if(cross == "DO") {
#    sanger = cbind("A_J" = geno(sanger)$GT[,1,drop = FALSE], 
#             "C57BL_6J" = "0/0",
#             geno(sanger)$GT[,2:7,drop = FALSE])
#  } else if(cross == "HS") {
    sanger = cbind(geno(sanger)$GT[,1:4,drop = FALSE],
             "C57BL_6J" = "0/0",
             geno(sanger)$GT[,5:7,drop = FALSE])
#  } # else

  # Convert allele calls to numeric values.
  sanger = (sanger != "0/0") * 1

  # Make the MAF between 1/8 and 4/8.
  flip = which(rowSums(sanger) > 4)
  sanger[flip,] = 1 - sanger[flip,,drop = FALSE]
  rm(flip)

###DMG
### I'm moving this outside of the function.
  # Create the survival object.
#  surv = Surv(pheno$days, pheno$cataract)

  # Null model.
###DMG
### Put the null logistic regression or linear model here.
  null.mod = coxph(surv ~ addcovar)
  null.ll = logLik(null.mod)
  pv = rep(0, nrow(sanger))

  # Get the unique SDPs between each pair of markers and
  # calculate the COXPH LOD.

  # CoxPH function.
  coxph.fxn = function(snp.rng, local.probs) {

    # Get the SDPs.
    sdp.nums = sanger[snp.rng,] %*% 2^(7:0)
    sdps2keep = which(!duplicated(sdp.nums))
    cur.sdps = sanger[snp.rng,,drop = FALSE][sdps2keep,,drop = FALSE]
    unique.sdp.nums = sdp.nums[sdps2keep]
    m = match(sdp.nums, unique.sdp.nums)

    # Multiply the SDPs by the haplotype probabilities.
    cur.alleles = tcrossprod(cur.sdps, local.probs)
    cur.ll = rep(null.ll, nrow(cur.sdps))

    # Check for low allele frequencies and remove SDPs with too
    # few samples carrying one allele.
    sdps.to.use = which(rowSums(cur.alleles) > 0.5)

    sex.col = which(colnames(addcovar) == "sex")
    if(length(sex.col) != 1) {
      stop("One of the columns of addcovar MUST be named 'sex'.")
    } # if(length(sex.col) != 1)

    # Run the Cox PH model at each unique SDP.
    for(j in sdps.to.use) {

###DMG
### Put the logistic regression or linear model here.

      # For the X chromosome we map with sex as an interactive
      # covariate with genotype.
      mod = coxph(surv ~ addcovar + cur.alleles[j,] + 
                  addcovar[,sex.col] * cur.alleles[j,])
      cur.ll[j] = logLik(mod)

    } # for(j)

    # This is the LRS.
    cur.ll = cur.ll - null.ll

    # Return the results.
    cur.ll[m]

   } # coxph.fxn()

   # SNPs before the first marker.
   snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
   if(length(snp.rng) > 0) {

     pv[snp.rng] = coxph.fxn(snp.rng, obj$probs[,,1])

   } # if(length(snp.rng) > 0)

   # SNPs between Markers.
   for(i in 1:(nrow(obj$markers)-1)) {

     snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                     sanger.hdr$POS <= obj$markers[i+1,3])

     if(length(snp.rng) > 0) {

       # Take the mean of the haplotype probs at the surrounding markers.
       pv[snp.rng] = coxph.fxn(snp.rng, (obj$probs[,,i] + 
                     obj$probs[,,i+1]) * 0.5)

     } # if(length(snp.rng) > 0)

   } # for(i)

  # SNPs after the last marker.
  snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
  if(length(snp.rng) > 0) {

    pv[snp.rng] = coxph.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])

  } # if(length(snp.rng) > 0)

  # Convert LRS to p-values using the chi-squared distribution.
  # Note that we have more degrees of freedom in the model.
  pv = pchisq(2 * pv, df = 13, lower.tail = FALSE)
  pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)

  save(pv, file = paste0(file.prefix, "_chr", chr, ".Rdata"))

  png(paste0(file.prefix, "_chr", chr,".png"), width = 2000, 
      height = 1600, res = 200)
  plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
  mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
  dev.off()

  # Return the positions and p-values.
  return(pv)

} # workfxn.xchr()


# Set up the worker cluster.
cl = makeCluster(ncl)
registerDoParallel(cl)
tmp = clusterEvalQ(cl, library(DOQTL))
tmp = clusterEvalQ(cl, library(VariantAnnotation))
tmp = clusterEvalQ(cl, library(regress))
tmp = clusterEvalQ(cl, library(survival))
### DMG
### Change this to export the "surv" object.
#clusterExport(cl, c("pheno", "addcovar", "snp.file", "outdir", "cross"))
clusterExport(cl, c("surv", "addcovar", "snp.file", "outdir", "cross"))

### DMG
### Comment out for now...
#result = foreach(i = iter(data)) %dopar% {

#  workfxn(i)

#} # for(each(i)

#save(result, file = paste0(file.prefix, ".Rdata"))

#stopCluster(cl)

##W/o cluster

result = vector("list", length(data))
names(result) = names(data)
for(i in 1:19) {
  print(i)
  result[[i]] = workfxn(data[[i]])
} #for(i)
print("X")
result[["X"]] = workfxn.xchr(data[["X"]])

save(result, file = paste0(file.prefix, ".Rdata"))

  
# Plotting function.
setwd(outdir)
files = dir(pattern = file.prefix)
files = files[files != paste0(file.prefix, ".Rdata")]
png.files = grep("png$", files)
if(length(png.files) > 0) {
  files = files[-png.files]
}
num = gsub(paste0("^", file.prefix, "_chr|\\.Rdata$"), "", files)
files = files[order(as.numeric(num))]

data = vector("list", length(files))
names(data) = num[order(as.numeric(num))]
for(i in 1:length(files)) {

  print(i)
  load(files[i])
  data[[i]] = pv
  data[[i]][,6] = -log10(data[[i]][,6])

} # for(i)

num.snps = sapply(data, nrow)
chrs = c(1:19, "X")

xlim = c(0, sum(num.snps))
ylim = c(0, max(sapply(data, function(z) { max(z[,6]) })))

# This plots all chromosomes.
chrlen = get.chr.lengths()[1:20]
chrsum = cumsum(chrlen)
chrmid = c(1, chrsum[-length(chrsum)]) + chrlen * 0.5
names(chrmid) = names(chrlen)

png(paste0(file.prefix, "_QTL.png"), width = 2000, height = 1600, res = 200)
plot(-1, -1, col = 0, xlim = c(0, max(chrsum)), ylim = ylim, xlab = "",
     ylab = "-log10(p-value)", las = 1, main = plot.title, xaxt = "n")
for(i in 1:length(data)) {
  print(i)
  pos = data[[i]][,3] * 1e-6 + c(0, chrsum)[i]
  points(pos, data[[i]][,6], col = c("black", "grey50")[i %% 2 + 1],
         pch = 20)
} # for(i)
mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.5)
dev.off()
