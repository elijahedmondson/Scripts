# Lifespan, females only, COXPH with unique SDPs.
library(DOQTL)
library(doParallel)
library(foreach)
library(Rsamtools)
library(VariantAnnotation)
library(GenomicRanges)
library(survival)
library(regress)
library(Rmpi)
options(stringsAsFactors = F)
setwd("/Users/elijah/Desktop/R/QTL/WD")

# Pass in the number of clusters (nodes).
args = commandArgs(trailingOnly = TRUE)
ncl = as.numeric(args[[1]])

####################
# THINGS TO CHANGE

# Set the output file directory.
outdir = "/Users/elijah/Desktop/R/QTL/WD/hq_snps"

# Load in your data. This file contains pheno, probs, markers and K.
HZE <- read.csv("~/Desktop/HZE-Fe.csv")
load(file = "/Users/elijah/Desktop/R/QTL/WD/model.probs.Rdata")
probs <- model.probs
rm(model.probs)

load(file = "/Users/elijah/Desktop/R/QTL/WD/K.Rdata")

pheno = data.frame(row.names = HZE$row.names, sex = as.numeric(HZE$sex == "M"),  
                   days = as.numeric(HZE$Cataract.2.0.Score),
                   cataract = as.numeric(HZE$Cataract.2.0.Event),
                   LSA = as.numeric(HZE$Lymphoma))
rm(HZE)
covar = data.frame(sex = as.numeric(pheno$sex == "M"))
addcovar = cbind(sex = as.numeric(factor(pheno$sex)) - 1)
rownames(addcovar) = rownames(pheno)
rownames(covar) = rownames(pheno)

load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
markers <- MM_snps
rm(MM_snps)
snps = markers[markers[,1] %in% dimnames(probs)[[3]],]
markers <- snps
rm(snps)

stopifnot(nrow(markers) == dim(probs)[3])
stopifnot(markers[,1] == dimnames(probs)[[3]])

# Set the Sanger SNP file location.
sdp.file = "/Users/elijah/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz"
snp.file = "ftp://ftp.jax.org/SNPtools/variants/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"

# Set the cross type (DO or HS).
cross = "HS"

# Set a file prefix for the output files.
file.prefix = "Cataract_Latency_coxph"

# Set the plot title.
plot.title = "Cataract Latency, CoxPH, HQ SNPs"

####################

# Set up covariates
rownames(addcovar) = rownames(pheno)

## Remove unused samples (?)
pheno = pheno[rownames(pheno) %in% dimnames(probs)[[1]],,drop = FALSE]
probs = probs[dimnames(probs)[[1]] %in% rownames(pheno),,]
probs = probs[match(rownames(pheno), dimnames(probs)[[1]]),,]
probs = probs[,,dimnames(probs)[[3]] %in% markers[,1]]

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

setwd(outdir)

# Make a function for each worker to execute.
workfxn = function(data) {

  chr = data$markers[1,2]

  setwd(outdir)

  # Get the Sanger SNPs.
  strains = sub("/", "_", do.colors[,2])
  if(cross == "HS") {
    strains = sub("/", "_", hs.colors[,2])
  } # if(cross = "HS")
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
  sanger.hdr = data.frame(ID = names(rowData(sanger)), CHR = as.character(seqnames(sanger)),
               POS = start(sanger), REF = as.character(fixed(sanger)$REF),
               ALT = alt, stringsAsFactors = FALSE)
  rm(alt)

  # Add C57BL/6J to the Sanger SNPs.
  if(cross == "DO") {
    sanger = cbind("A_J" = geno(sanger)$GT[,1,drop = FALSE], 
             "C57BL_6J" = "0/0",
             geno(sanger)$GT[,2:7,drop = FALSE])
  } else if(cross == "HS") {
    sanger = cbind(geno(sanger)$GT[,1:4,drop = FALSE],
             "C57BL_6J" = "0/0",
             geno(sanger)$GT[,5:7,drop = FALSE])
  } # else

  # Convert allele calls to numeric values.
  sanger = (sanger != "0/0") * 1

  # Make the MAF between 1/8 and 4/8.
  flip = which(rowSums(sanger) > 4)
  sanger[flip,] = 1 - sanger[flip,,drop = FALSE]
  rm(flip)

  # Create the survival object.
  surv = Surv(pheno$days, pheno$cataract)

  # Null model.
  null.mod = coxph(surv ~ addcovar)
  null.ll = logLik(null.mod)
  pv = rep(0, nrow(sanger))

  # Get the unique SDPs between each pair of markers and
  # calculate the COXPH LOD.

  # CoxPH function.
  coxph.fxn = function(snp.rng, local.probs) {

    # Get the SDPs.
    sdps = sanger[snp.rng,] %*% 2^(7:0)
    sdps2keep = which(!duplicated(sdps))
    cur.sdps = sanger[snp.rng,,drop = FALSE][sdps2keep,,drop = FALSE]
    unique.sdps = sdps[sdps2keep]
    m = match(sdps, unique.sdps)

    # Calculate the SNP probs.
    cur.alleles = tcrossprod(cur.sdps, local.probs)
    cur.ll = rep(0, nrow(cur.sdps))

    # Run the Cox PH model at each unique SDP.
    for(j in 1:nrow(cur.sdps)) {

      mod = coxph(surv ~ addcovar + cur.alleles[j,])
      cur.ll[j] = logLik(mod)

    } # for(j)
    # This is the LRS.
    cur.ll = cur.ll - null.ll

    # Return the results.
    cur.ll[m]

   } # coxph.fxn()

   # SNPs before the first marker.
   snp.rng = which(sanger.hdr$POS <= data$markers[1,3])
   if(length(snp.rng) > 0) {

     pv[snp.rng] = coxph.fxn(snp.rng, data$probs[,,1])

   } # if(length(snp.rng) > 0)

   # SNPs between Markers.
   for(i in 1:(nrow(data$markers)-1)) {

     snp.rng = which(sanger.hdr$POS > data$markers[i,3] &
                     sanger.hdr$POS <= data$markers[i+1,3])

     if(length(snp.rng) > 0) {

       # Take the mean of the haplotype probs at the surrounding markers.
       pv[snp.rng] = coxph.fxn(snp.rng, (data$probs[,,i] + 
                     data$probs[,,i+1]) * 0.5)

     } # if(length(snp.rng) > 0)

   } # for(i)

  # SNPs after the last marker.
  snp.rng = which(sanger.hdr$POS > data$markers[nrow(data$markers),3])
  if(length(snp.rng) > 0) {

    pv[snp.rng] = coxph.fxn(snp.rng, data$probs[,,nrow(data$markers)])

  } # if(length(snp.rng) > 0)

  # Convert LRS to p-values using the chi-squared distribution.
  pv = pchisq(2 * pv, df = 1, lower.tail = FALSE)
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

# Set up the worker cluster.
cl = makeCluster(ncl, type = "MPI")
registerDoParallel(cl)
tmp = clusterEvalQ(cl, library(DOQTL))
tmp = clusterEvalQ(cl, library(VariantAnnotation))
tmp = clusterEvalQ(cl, library(regress))
tmp = clusterEvalQ(cl, library(survival))
clusterExport(cl, c("pheno", "addcovar", "snp.file", "outdir", "cross"))

result = foreach(i = iter(data)) %dopar% {

  workfxn(i)

} # for(each(i)

save(result, file = paste0(file.prefix, ".Rdata"))

stopCluster(cl)

################# Without cluster
result1 <- workfxn(data[[1]])
result2 <- workfxn(data[[2]])
result3 <- workfxn(data[[3]])
result4 <- workfxn(data[[4]])
result5 <- workfxn(data[[5]])
result6 <- workfxn(data[[6]])
result7 <- workfxn(data[[7]])
result8 <- workfxn(data[[8]])
result9 <- workfxn(data[[9]])
result10 <- workfxn(data[[10]])
result11 <- workfxn(data[[11]])
result12 <- workfxn(data[[12]])
result13 <- workfxn(data[[13]])
result14 <- workfxn(data[[14]])
result15 <- workfxn(data[[15]])
result16 <- workfxn(data[[16]])
result17 <- workfxn(data[[17]])
result18 <- workfxn(data[[18]])
result19 <- workfxn(data[[19]])
resultX <- workfxn(data[["X"]])


result <- for(i=1:length(data)) {
  workfxn(data[[i]])
}
  
  
# Plotting function.
setwd(outdir)
files = dir(pattern = file.prefix)
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
png(paste0(file.prefix, "_QTL.png"), width = 2000, height = 1600, res = 200)
plot(-1, -1, col = 0, xlim = xlim, ylim = ylim, xlab = "",
     ylab = "-log10(p-value)", las = 1, main = plot.title)
idx = 1
for(i in 1:length(data)) {
  print(i)
  rng = idx:(idx + num.snps[i] - 1)
  points(rng, data[[i]][,6], col = c("black", "grey50")[i %% 2 + 1], pch = 20)
  idx = idx + num.snps[i]
} # for(i)
dev.off()
