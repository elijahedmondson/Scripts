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

#### Files ####
####
load(file = "C:/Users/edmondsonef/Desktop/QTL/HZEproject.Rdata")
load(file = "C:/Users/edmondsonef/Desktop/QTL/HS.colors.Rdata")
sdp.file = "C:/Users/edmondsonef/Desktop/QTL/HS_Sanger_SDPs.txt.bgz"
#load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
outdir = "C:/Users/edmondsonef/Desktop/R-plots/"


#### Phenotype ####
####

total <- read_excel("C:/Users/edmondsonef/Desktop/Cataract/CATARACT_final.xlsx", sheet ="simple")
pheno <- data.frame(row.names = total$row.names,
                    sex = as.numeric(total$sex == "M"),
                    group = as.character(total$groups),
                    albino = as.numeric(total$`coat color` == "albino"),
                    cat_average = as.numeric(total$`cat_average`),
                    cat_sum = as.numeric(total$`cat_sum`),
                    cat_difference = as.numeric(total$`cat_difference`))
HZE <- subset(pheno, group == "HZE")
Gamma <- subset(pheno, group == "Gamma")
Un <- subset(pheno, group == "Unirradiated")
Allirr <- subset(pheno, group != "Unirradiated")

addcovar <- matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))




#### ASSOCIATION MAPPING ####
####
####
# qtl <- GRSD.assoc(pheno = HZE, pheno.col = 4, probs, K, addcovar, markers, snp.file,
#                   outdir = "C:/Users/edmondsonef/Desktop/R-plots/", tx = "HZE",
#                   sanger.dir = "C:/Users/edmondsonef/Desktop/QTL/HS.sanger.files/")

qtl <- GRSD.gauss(pheno = HZE, pheno.col = 4, probs, K, addcovar, markers, snp.file,
                  outdir = "C:/Users/edmondsonef/Desktop/R-plots/", tx = "HZE",
                  sanger.dir = "C:/Users/edmondsonef/Desktop/QTL/HS.sanger.files/")

# qtl <- GRSD.poisson(pheno = HZE, pheno.col = 4, probs, K, addcovar, markers, snp.file,
#                     outdir = "C:/Users/edmondsonef/Desktop/R-plots/", tx = "HZE",
#                     sanger.dir = "C:/Users/edmondsonef/Desktop/QTL/HS.sanger.files/")

# qtl <- scanone(pheno = Gamma, pheno.col = 4, probs = probs, K = K, addcovar = addcovar, snps = markers)
# 
# 
# qtl <- scanone.assoc(pheno = HZE, pheno.col = 4, probs = probs, K = K, cross = HS,
#                      addcovar = addcovar, markers = markers, sdp.file = sdp.file, ncl = 4)



#### FIND THE MAX LOD ####

max.LOD = result[[7]]$ID[which(-log10(result[[7]]$pv) > 10)]
max.LOD

max.LOD.position <- result[[7]]$POS[which(-log10(result[[7]]$pv) > 100)]
max.LOD.position

mgi = get.mgi.features(chr = 7, start = 87491804, end = 88689056, type = "gene", source = "MGI")




#### Permutation analysis & Significatn thresholds ####
####
# perms <- Scanone.assoc.perms(perms = 200, pheno = Gamma, pheno.col = "AML.t", probs = model.probs, 
#                              K = K, tx = "Gamma", addcovar = addcovar, markers = MM_snps, 
#                              sdp.file = sdp.file, ncl = 4)

perms <- GRSDgauss.perms(perms = 200, chr = 1:19, Xchr = FALSE,
                         pheno = HZE, pheno.col = 4, probs, K, addcovar,
                         markers, snp.file, outdir = "C:/Users/edmondsonef/Desktop/R-plots/",
                         tx = "", sanger.dir = "C:/Users/edmondsonef/Desktop/QTL/HS.sanger.files/")

get.sig.thr(perms, alpha = 0.01, Xchr = F)
thr1 = quantile(perms, probs = 0.90)
thr2 = quantile(perms, probs = 0.95)
thr3 = quantile(perms, probs = 0.99)




#### Plotting ####
####

plot.scanone.assoc(res, bin.size = 100, main = "")

png("HZE_scanoneassoc_cat_sum.png", width = 2400, height = 1080, res = 200)
plot.scanone.assoc(res, bin.size = 100)
dev.off()

png("_CHR_4.png", width = 2000, height = 1600, res = 128)
plot.scanone.assoc(qtl, chr = 7, bin.size = 100)
dev.off()

layout(matrix(3:1, 3, 1))
par(mfrow = c(2,2), mar=c(1, 4, 1, 1) + 0.1)
loop.hs.qtl(HZE.cat2, chr=17, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
loop.hs.qtl(Gamma.cat2, chr=17, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
loop.hs.qtl(Unirradiated.cat2, chr=17, bin.size = 100, main = "Unirradiated", ylim=c(0,15))
loop.hs.qtl(Allirr.cat2, chr=17, bin.size = 100, main = "All irradiated", ylim=c(0,15))


# layout(matrix(3:1, 3, 1))
# par(mfrow = c(2,2), mar=c(1, 4, 1, 1) + 0.1)
# loop.hs.qtl(HZE.cat2, chr=17, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
# DOQTL:::plot.scanone.assoc(Gamma.cat2, chr=17, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
# DOQTL:::plot.scanone.assoc(Unirradiated.cat2, chr=17, bin.size = 100, main = "Unirradiated", ylim=c(0,15))
# DOQTL:::plot.scanone.assoc(Allirr.cat2, chr=17, bin.size = 100, main = "All irradiated", ylim=c(0,15))
# 
# par(mfrow = c(3,1), mar=c(1, 4, 1, 1) + 0.5)
# DOQTL:::plot.scanone.assoc(HZE.days, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
# abline(a = 13, b = 0, col = "red")
# DOQTL:::plot.scanone.assoc(Gamma.days, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
# abline(a = 13, b = 0, col = "red")
# DOQTL:::plot.scanone.assoc(Unirradiated.days, bin.size = 100, main = "Unirradiated", ylim=c(0,15))
# abline(a = 13, b = 0, col = "red")




#### LINKAGE MAPPING ####
####
# perms = scanone.perm(pheno = pheno, pheno.col = "AML", probs = model.probs, addcovar = addcovar,
#                      snps = MM_snps, path = "~/Desktop/",
#                      nperm = 100)
# get.sig.thr(perms, alpha = 0.01, Xchr = TRUE)
# thr1 = quantile(perms, probs = 0.90)
# thr2 = quantile(perms, probs = 0.95)
# thr3 = quantile(perms, probs = 0.99)
# 
# plot(qtl, sig.thr = c(thr1, thr2, thr3), main = "PSC")
# 


#### MGI features
####
# interval = bayesint(qtl, chr = 11)
# interval
# mgi = get.mgi.features(chr = interval[1,2], start = interval[1,3],
#                        end = interval[3,3], type = "gene", source = "MGI")
# mgi = get.mgi.features(chr = 14, start = 68772000, 
#                        end = 68774069, type = "gene", source = "MGI")
# nrow(mgi)
# head(mgi)
# 
# ma = assoc.map(pheno = pheno, pheno.col = "albino", probs = probs, K = K, addcovar = addcovar,
#                snps = MM_snps, chr = interval[1,2], start = interval[1,3], end = interval[3,3])
# coefplot(qtl, chr = 7, cross = "HS", colors = "HS")
# tmp = assoc.plot(ma, thr = 1)
# unique(tmp$sdps)
# 



#### Helper Functions ####
####
################################################################################
loop.hs.qtl = function(qtl, title, bin.width = 1000, ...) {
  library(GenomicRanges)
  library(BSgenome.Mmusculus.UCSC.mm10)
  new.qtl = NULL
  for(chr in 1:length(qtl)) {
    
    print(chr)
    
    # Create SNP bins with given bin.width
    brks = cut(x = 1:length(qtl[[chr]]), breaks = length(qtl[[chr]]) / bin.width)
    # Split up the SNP positions and get the mean.
    pos = split(start(qtl[[chr]]), brks)
    pos = sapply(pos, mean)
    # Split up the p-values and get the max.
    pv = split(mcols(qtl[[chr]])$p.value, brks)
    pv = sapply(pv, min)
    
    # Make a single new GRanges object to return.
    gr = GRanges(seqnames = seqnames(qtl[[chr]])[1],
                 ranges = IRanges(start = pos, width = 1), p.value = pv)
    
    if(chr == 1) {
      new.qtl = gr
    } else {
      new.qtl = c(new.qtl, gr)
    } # else
    
  } # for(chr)
  
  # Get the chromosome lengths.
  chrlen = seqlengths(BSgenome.Mmusculus.UCSC.mm10)
  names(chrlen) = sub("^chr", "", names(chrlen))
  chrlen = chrlen[seqlevels(new.qtl)] * 1e-6
  
  # Add the chr lengths to the chromosomes for plotting.
  # Switch positions to genome Mb.
  gmb = start(new.qtl) * 1e-6
  for(chr in 2:length(chrlen)) {
    
    wh = which(seqnames(new.qtl) == names(chrlen)[chr])
    gmb[wh] = gmb[wh] + sum(chrlen[1:(chr - 1)])
    
  } # for(chr)
  
  # Get chromosome mid-points for plotting the Chr name.
  chrmid = (chrlen / 2) + cumsum(c(1, chrlen[-length(chrlen)]))
  
  # Make the plot.
  col = rep(rgb(0,0,0), length(new.qtl))
  even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
  col[even.chr] = rgb(0.7,0.7,0.7)
  plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
       col = col, las = 1, xlab = "", ylab = "-log10(p-value)", main = title)
  mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)
  
  return(new.qtl)
  
}
################################################################################
plot.scanone.assoc = function(x, chr, bin.size = 1000, sig.thr,
                              sig.col = "red", ...) {
  
  if(!missing(chr)) {
    x = x[names(x) %in% chr]
  } # if(!missing(chr))
  
  chrlen = get.chr.lengths()
  chrlen = chrlen[names(chrlen) %in% names(x)]
  chrsum = cumsum(chrlen)
  chrmid = c(0, chrsum[-length(chrsum)]) + diff(c(0, chrsum)) * 0.5
  names(chrmid) = names(chrsum)
  
  if(missing(chr)) {
    autosomes = names(x)[which(!is.na(as.numeric(names(x))))]
    chr = factor(names(x), levels = c(autosomes, "X", "Y", "M"))
  } # if(!missing(chr))
  
  pos = lapply(x, start)
  pv =  lapply(x, mcols)
  pv =  lapply(pv, "[[", 1)
  
  for(c in 1:length(x)) {
    
    bins = seq(1, length(pv[[c]]), bin.size)
    bins[length(bins) + 1] = length(pv[[c]])
    pos2 = rep(0, length(bins))
    pv2  = rep(0, length(bins))
    
    for(i in 1:(length(bins)-1)) {
      wh = which.min(pv[[c]][bins[i]:bins[i+1]])
      wh = wh + bins[i] - 1
      pos2[i] = pos[[c]][wh] * 1e-6 + max(0, chrsum[c - 1])
      pv2[i]  = pv[[c]][wh]
    } # for(i)
    
    pos[[c]] = pos2
    pv[[c]]  = -log(pv2, 10)
    
  } # for(c)
  
  # If we are plotting more than one chormosome, color alternate
  # chromosomes grey and black.
  col = 1
  chr = rep(chr, sapply(pos, length))
  if(length(pos) > 1) {
    col = as.numeric(chr) %% 2 + 1
  } # if(length(chr) > 1)
  
  plot(unlist(pos), unlist(pv), pch = 16, col = c("black", "grey60")[col],
       las = 1, xaxt = "n", xlab = "", ylab = "-log10(p-value)", xaxs = "i", ...)
  mtext(text = names(chrmid), side = 1, line = 2.5, at = chrmid, cex = 2)
  
  if(length(pos) == 1) {
    
    axis(side = 1)
    
  } # if(length(pos) == 1)
  
  if(!missing(sig.thr)) {
    
    add.sig.thr(sig.thr = sig.thr, sig.col = sig.col, chrsum = chrsum)
    
  } # if(!missing(sig.thr))
  
} # plot.scanone.assoc()
################################################################################
scanone.assoc = function(pheno, pheno.col, probs, K, addcovar, intcovar, markers,
                         cross = c("DO", "CC", "HS"), sdp.file, ncl) {
  cl = makeCluster(ncl)
  registerDoParallel(cl)
  # Synch up markers and haplotype probs.
  markers = markers[!is.na(markers[,3]),]
  markers = markers[markers[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,dimnames(probs)[[3]] %in% markers[,1]]
  # Put the marker positions on a Mb scale.
  if(any(markers[,3] > 200, na.rm = TRUE)) {
    markers[,3] = markers[,3] * 1e-6
  } # if(any(markers[,3] > 200)
  # Synch up the sample IDs.
  tmp = synch.sample.IDs(pheno = pheno, probs = probs, K = K, addcovar = addcovar)
  pheno = tmp$pheno
  addcovar = tmp$addcovar
  probs = tmp$probs
  K = tmp$K
  rm(tmp)
  gc()
  # Get the unique chromosomes.
  chr = unique(markers[,2])
  chr = lapply(chr, "==", markers[,2])
  chr = lapply(chr, which)
  names(chr) = unique(markers[,2])
  # Split up the data and create a list with elements for each 
  # chromosome.
  data = vector("list", length(chr))
  for(i in 1:length(chr)) {
    data[[i]] = list(pheno = pheno, pheno.col = pheno.col, 
                     addcovar = addcovar, probs = probs[,,chr[[i]]], K = K[[i]], 
                     markers = markers[chr[[i]],])
  } # for(i)
  names(data) = names(chr)
  rm(probs, markers, K, pheno, pheno.col, addcovar)
  # Load the required libraries on the cores.
  clusterEvalQ(cl = cl, expr = library(DOQTL))
  clusterEvalQ(cl = cl, expr = library(Rsamtools))
  clusterEvalQ(cl = cl, expr = library(regress))
  res = foreach(obj = iter(data)) %dopar% {
    
    s1.assoc(obj, sdp.file)
  } # foreach(c)
  names(res) = names(data)
  stopCluster(cl)
  class(res) = c("scanone.assoc", class(res))
  return(res)
} # scanone.assoc()
# Helper function to perform association mapping on one autosome.
# obj: data object of the type created in scanone.assoc().
# sdp.file: Tabix file created by condense.sanger.snps().
################################################################################
s1.assoc = function(obj, sdp.file) {
  # Get the samples that are not NA for the current phenotype.
  sample.keep = which(!is.na(obj$pheno[,obj$pheno.col]) & 
                        rowSums(is.na(obj$addcovar)) == 0)
  # Calculate variance components and change each kinship matrix to be the
  # error covariance correction matrix.
  mod = regress(obj$pheno[,obj$pheno.col] ~ obj$addcovar, ~obj$K,
                pos = c(TRUE, TRUE))
  obj$K = mod$sigma[1] * obj$K + mod$sigma[2] * diag(nrow(obj$K))
  rm(mod)
  # Read in the unique SDPs.
  tf = TabixFile(sdp.file)
  sdps = scanTabix(file = sdp.file, param = GRanges(seqnames = obj$markers[1,2],
                                                    ranges = IRanges(start = 0, end = 200e6)))[[1]]
  sdps = strsplit(sdps, split = "\t")
  sdps = matrix(unlist(sdps), ncol = 3, byrow = T)
  chr  = sdps[1,1]
  pos  = as.numeric(sdps[,2])
  sdps = as.numeric(sdps[,3])
  # Create a matrix of SDPs.
  sdp.mat = matrix(as.numeric(intToBits(1:2^8)), nrow = 32)
  sdp.mat = sdp.mat[8:1,]
  dimnames(sdp.mat) = list(LETTERS[1:8], 1:2^8)
  # Between each pair of markers, get the unique SDPs and their genomic
  # positions. Use the DO genoprobs to create DO genotypes.
  # Get a set of overlaps between the markers and SDP positions.
  sdp.gr = GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1))
  # Include 0 and 200 Mb to capture SNPs before the first markers and 
  # after the last marker.
  markers.gr = GRanges(seqnames = obj$markers[1,2], ranges = IRanges(
    start = c(0, obj$markers[,3]) * 1e6, 
    end = c(obj$markers[,3], 200) * 1e6))
  ol = findOverlaps(query = markers.gr, subject = sdp.gr)
  ol = split(subjectHits(ol), queryHits(ol))
  probs.idx = as.numeric(names(ol))
  unique.sdps = lapply(ol, function(z) { unique(sdps[z]) })
  num.sdps = sum(sapply(unique.sdps, length))
  geno = matrix(0, nrow = nrow(obj$pheno), ncol = num.sdps,
                dimnames = list(rownames(obj$pheno), 1:num.sdps))
  # This maps the positions in the geno matrix back to the genomic positions
  # of the SDPs.
  map = rep(0, length(pos))
  idx = 0
  # Start of chromosome, sdps before the first marker.
  if(probs.idx[1] == 1) {
    i = 1
    # Get the range of SDPs.
    rng = (idx + 1):(idx + length(unique.sdps[[i]]))
    # Multiply the first genoprobs by the SDPs before the first marker.
    geno[,rng] = obj$probs[,,probs.idx[i]] %*% sdp.mat[,unique.sdps[[i]]]
    # Place the SDPs in the map.
    map[ol[[i]]] = match(sdps[ol[[i]]], unique.sdps[[i]]) + idx
    # Increment the index.
    idx = idx + length(unique.sdps[[i]])
  } # if(probs.idx[1] == 1) 
  # SDPs bracketed by two markers.
  wh = which(probs.idx > 1 & probs.idx <= dim(obj$probs)[3])
  for(i in wh) {
    rng = (idx + 1):(idx + length(unique.sdps[[i]]))
    # Use the mean genoprobs between two markers and multiply by the SDPs.
    geno[,rng] = 0.5 * (obj$probs[,,probs.idx[i] - 1] + obj$probs[,,probs.idx[i]]) %*% 
      sdp.mat[,unique.sdps[[i]]]
    map[ol[[i]]] = match(sdps[ol[[i]]], unique.sdps[[i]]) + idx
    idx = idx + length(unique.sdps[[i]])
  } # for(i)
  # End of chromosome, sdps after the last marker.
  i = length(probs.idx)
  if(probs.idx[i] > dim(obj$probs)[3]) {
    rng = (idx + 1):(idx + length(unique.sdps[[i]]))
    geno[,rng] = obj$probs[,,dim(obj$probs)[3]] %*% sdp.mat[,unique.sdps[[i]]]
    map[ol[[i]]] = match(sdps[ol[[i]]], unique.sdps[[i]]) + idx
  } # if(probs.idx[length(probs.idx)] > dim(obj$probs)[3])
  r2 = 0
  # X Chromosome, separate females and males and then combine.
  if(chr == "X") {
    # Verify that sex is one of the covariates.
    if(!any("addcovar" == names(obj))) {
      stop(paste("In order to map on the X chromosome, you must supply the sex",
                 "of each sample in \'addcovar\', even if they are all the same sex."))
    } # if(!any("addcovar" == names(obj)))
    if(length(grep("sex", colnames(obj$addcovar), ignore.case = T)) == 0) {
      stop(paste("In order to map on the X chromosome, you must supply the sex",
                 "of each sample in \'addcovar\', even if they are all the same sex."))
    } # if(grep("sex", colnames(obj$addcovar), ignore.case = T))
    # Get the sex of each sample.
    sex.col = grep("sex", colnames(obj$addcovar), ignore.case = TRUE)
    females = which(obj$addcovar[,sex.col] == 0)
    males   = which(obj$addcovar[,sex.col] == 1)
    if(length(females) == 0 & length(males) == 0) {
      stop(paste("Sex is not coded using 0 for female and 1 for males. Please",
                 "set the sex column in addcovar to 0 for females and 1 for males."))
    } # if(length(females) == 0 & length(males) == 0)
    # Separate the male and female genotypes in the same matrix. This will be 
    # a block matrix with all zeros for females in rows with male samples and
    # all zeros for males in rows with female samples.
    newgeno = matrix(0, nrow(geno), 2 * ncol(geno))
    newgeno[females,1:ncol(geno)] = geno[females,]
    newgeno[males,(ncol(geno)+1):ncol(newgeno)] = geno[males,]
    r2 = matrixeqtl.snps(pheno = obj$pheno[sample.keep,obj$pheno.col,drop = FALSE],
                         geno = newgeno[sample.keep,,drop = FALSE],
                         K = obj$K[sample.keep, sample.keep,drop = FALSE],
                         addcovar = obj$addcovar[sample.keep,,drop = FALSE])
    r2 = r2[1:ncol(geno)] + r2[(ncol(geno)+1):length(r2)]
    rm(newgeno)
  } else {
    # Autosomes.
    # Calculate the R^2 for each SDP.
    r2 = matrixeqtl.snps(pheno = obj$pheno[sample.keep,obj$pheno.col,drop = FALSE],
                         geno = geno[sample.keep,,drop = FALSE],
                         K = obj$K[sample.keep,sample.keep,drop = FALSE],
                         addcovar = obj$addcovar[sample.keep,,drop = FALSE])
  } # else
  # Convert R^2 to LRS.
  lrs = -length(sample.keep) * log(1.0 - r2)
  # Convert the LRS to p-values.
  pv = pchisq(lrs, df = 1, lower.tail = FALSE)
  # Place the results in the correct locations and return.
  return(GRanges(seqnames = Rle(chr, length(pos)), ranges = IRanges(start = pos, 
                                                                    width = 1), p.value = pv[map]))
} # s1.assoc()
################################################################################
# Plotting for scanone.assoc.
plot.scanone.assoc = function(x, chr, bin.size = 1000, sig.thr, 
                              sig.col = "red", show.chr = TRUE, ...)  {
  if(!missing(chr)) {
    x = x[names(x) %in% chr]
  } # if(!missing(chr))
  chrlen = get.chr.lengths()
  chrlen = chrlen[names(chrlen) %in% names(x)]
  chrsum = cumsum(chrlen)
  chrmid = c(0, chrsum[-length(chrsum)]) + diff(c(0, chrsum)) * 0.5
  names(chrmid) = names(chrsum)
  if(missing(chr)) {
    autosomes = names(x)[which(!is.na(as.numeric(names(x))))]
    chr = factor(names(x), levels = c(autosomes, "X", "Y", "M"))
  } # if(!missing(chr))
  pos = lapply(x, start)
  pv =  lapply(x, mcols)
  pv =  lapply(pv, "[[", 1)
  for(c in 1:length(x)) {
    bins = seq(1, length(pv[[c]]), bin.size)
    bins[length(bins) + 1] = length(pv[[c]])
    pos2 = rep(0, length(bins))
    pv2  = rep(0, length(bins))
    for(i in 1:(length(bins)-1)) {
      wh = which.min(pv[[c]][bins[i]:bins[i+1]])
      wh = wh + bins[i] - 1
      pos2[i] = pos[[c]][wh] * 1e-6 + max(0, chrsum[c - 1])
      pv2[i]  = pv[[c]][wh]
    } # for(i)
    pos[[c]] = pos2
    pv[[c]]  = -log(pv2, 10)
  } # for(c)
  # If we are plotting more than one chormosome, color alternate 
  # chromosomes grey and black.
  col = 1
  chr = rep(chr, sapply(pos, length))
  if(length(pos) > 1) {
    col = as.numeric(chr) %% 2 + 1
  } # if(length(chr) > 1)
  plot(unlist(pos), unlist(pv), pch = 16, col = c("black", "grey60")[col],
       las = 1, xaxt = "n", xlab = "", ylab = "-log10(p-value)", xaxs = "i", ...)
  if(show.chr) {
    if(length(pos) == 1) {
      axis(side = 1)
      mtext(text = names(chrmid), side = 1, line = 2.5, at = chrmid, cex = 1.5)
    } else {
      mtext(text = names(chrmid), side = 1, line = 0.5, at = chrmid, cex = 1.5)
    } # else
  } # else
  if(!missing(sig.thr)) {
    add.sig.thr(sig.thr = sig.thr, sig.col = sig.col, chrsum = chrsum)
  } # if(!missing(sig.thr))
} # plot.scanone.assoc()
################################################################################
GRSD.assoc = function(pheno, pheno.col, probs, K, addcovar, markers, snp.file,
                      outdir = "~/Desktop/", tx = c("Gamma", "HZE", "Unirradiated", "All"),
                      sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/"){
  begin <- Sys.time()
  begin
  # COVARIATES #
  
  samples = intersect(rownames(pheno), rownames(probs))
  samples = intersect(samples, rownames(addcovar))
  samples = intersect(samples, rownames(K[[1]]))
  stopifnot(length(samples) > 0)
  print(paste("A total of", length(samples), tx, "samples are complete."))
  
  pheno = pheno[samples,,drop = FALSE]
  addcovar = addcovar[samples,,drop = FALSE]
  probs = probs[samples,,,drop = FALSE]
  
  
  
  # DEFINE TRAIT #
  
  file.prefix = paste(tx, pheno.col, sep = "_")
  
  plot.title = paste(tx, pheno.col, sep = " ")
  print(plot.title)
  
  trait = pheno[,pheno.col]
  print(table(trait))
  print(paste(round(100*(sum(trait) / length(samples)), digits = 1),
              "% display the", pheno.col, "phenotype in the", tx, "group."))
  
  # LOGISTIC REGRESSION MODEL #
  
  for(i in 1:length(K)) {
    K[[i]] = K[[i]][samples, samples]
  } # for(i)
  
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
  
  setwd(outdir)
  
  # MAPPING ANALYSES #
  
  result = vector("list", length(data))
  names(result) = names(data)
  print(paste("Mapping with", length(samples), tx, "samples..."))
  sanger.dir = sanger.dir
  
  for(i in 1:19) {
    print(paste("CHROMOSOME", i))
    timechr <- Sys.time()
    result[[i]] = GRSDbinom.fast(data[[i]], pheno, pheno.col, addcovar, tx, sanger.dir)
    print(paste(round(difftime(Sys.time(), timechr, units = 'mins'), digits = 2),
                "minutes..."))
  } #for(i)
  
  print("X CHROMOSOME")
  result[["X"]] = GRSDbinom.xchr.fast(data[["X"]], pheno, pheno.col, addcovar, tx, sanger.dir)
  
  print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
              "hours elapsed during mapping."))
  
  # Convert to GRangesList for storage
  chrs = c(1:19, "X")
  qtl = GRangesList(GRanges("list", length(result)))
  
  for(i in 1:length(chrs)) {
    print(i)
    qtl[[i]] <- GRanges(seqnames = Rle(result[[i]]$CHR),
                        ranges = IRanges(start = result[[i]]$POS, width = 1),
                        p.value = result[[i]]$pv)
  } # for(i)
  
  
  
  # PLOTTING
  plotter <- Sys.time()
  
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
  print("Plotting...")
  for(i in 1:length(files)) {
    
    load(files[i])
    data[[i]] = pv
    data[[i]][,6] = -log10(data[[i]][,6])
    
  } # for(i)
  
  num.snps = sapply(data, nrow)
  chrs = c(1:19, "X")
  
  xlim = c(0, sum(num.snps))
  ylim = c(0, max(sapply(data, function(z) { max(z[,6]) })))
  
  # PLOT ALL CHROMOSOMES #
  setwd(outdir)
  chrlen = get.chr.lengths()[1:20]
  chrsum = cumsum(chrlen)
  chrmid = c(1, chrsum[-length(chrsum)]) + chrlen * 0.5
  names(chrmid) = names(chrlen)
  
  png(paste0(file.prefix, "_QTL.png"), width = 2600, height = 1200, res = 200)
  plot(-1, -1, col = 0, xlim = c(0, max(chrsum)), ylim = ylim, xlab = "",
       ylab = "-log10(p-value)", las = 1, main = plot.title, xaxt = "n")
  for(i in 1:length(data)) {
    print(paste("Plotting chromosome", i))
    pos = data[[i]][,3] * 1e-6 + c(0, chrsum)[i]
    points(pos, data[[i]][,6], col = c("black", "grey50")[i %% 2 + 1],
           pch = 20)
  } # for(i)
  mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.5)
  dev.off()
  
  # Convert to GRangesList for storage
  chrs = c(1:19, "X")
  qtl = GRangesList(GRanges("list", length(result)))
  
  for(i in 1:length(chrs)) {
    print(i)
    qtl[[i]] <- GRanges(seqnames = Rle(result[[i]]$CHR),
                        ranges = IRanges(start = result[[i]]$POS, width = 1),
                        p.value = result[[i]]$pv)
  } # for(i)
  
  save(qtl, file.prefix, file = paste0(file.prefix, "_QTL.Rdata"))
  
  print(paste(round(difftime(Sys.time(), plotter, units = 'hours'), digits = 2),
              "hours elapsed during plotting."))
  
}
################################################################################
GRSDbinom.fast = function(obj, pheno, pheno.col, addcovar, tx, sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/") {
  chr = obj$markers[1,2]
  
  setwd(outdir)
  
  file.prefix = paste(tx, pheno.col, sep = "_")
  
  plot.title = paste(tx, pheno.col, sep = " ")
  
  strains = sub("/", "_", hs.colors[,2])
  
  load(file = paste0(sanger.dir, chr, ".Rdata"))
  
  null.mod = glm(pheno[,pheno.col] ~ addcovar, family = binomial(logit))
  #null.mod = glm(trait ~ addcovar, family = poisson(link = "log"))
  null.ll = logLik(null.mod)
  pv = rep(0, nrow(sanger))
  
  glm.fxn = function(snp.rng, local.probs) {
    
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
    sdps.to.use = which(rowSums(cur.alleles) > 1.0)
    
    # Run the model at each unique SDP.
    for(j in sdps.to.use) {
      
      
      # library(regress)
      # mod = regress(obj$pheno[,obj$pheno.col] ~ obj$addcovar, ~obj$K, pos = c(TRUE, TRUE))
      # obj$K = mod$sigma[1] * obj$K + mod$sigma[2] * diag(nrow(obj$K))
      # rm(mod)
      
      
      
      full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = binomial(logit))
      #full.mod = glm(trait ~ addcovar + cur.alleles[j,], family = poisson(link = "log"))
      cur.ll[j] = logLik(full.mod)
      
    } # for(j)
    
    # This is the LRS.
    cur.ll = cur.ll - null.ll
    
    # Return the results.
    cur.ll[m]
    
  } # glm.fxn()
  
  # SNPs before the first marker.
  snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,1])
    
  } # if(length(snp.rng) > 0)
  
  # SNPs between Markers.
  for(i in 1:(nrow(obj$markers)-1)) {
    
    snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                      sanger.hdr$POS <= obj$markers[i+1,3])
    
    if(length(snp.rng) > 0) {
      
      # Take the mean of the haplotype probs at the surrounding markers.
      pv[snp.rng] = glm.fxn(snp.rng, (obj$probs[,,i] +
                                        obj$probs[,,i+1]) * 0.5)
      
    } # if(length(snp.rng) > 0)
    
  } # for(i)
  
  # SNPs after the last marker.
  snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])
    
  } # if(length(snp.rng) > 0)
  
  # Convert LRS to p-values using the chi-squared distribution.
  pv = pchisq(2 * pv, df = 1, lower.tail = FALSE)
  pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)
  
  maxPV = max(-log10(pv[,6]))
  print(paste0("Maximum LOD: ", maxPV), digits = 1)
  rm(maxPV)
  
  save(pv, file = paste0(file.prefix, "_chr", chr, ".Rdata"))
  
  png(paste0(file.prefix, "_chr", chr,".png"), width = 2600,
      height = 1200, res = 130)
  plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
  mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
  dev.off()
  
  # Return the positions and p-values.
  return(pv)
  
  
} # GRSDbinom.fast()
################################################################################
GRSDbinom.xchr.fast = function(obj, pheno, pheno.col, addcovar, tx, sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/") {
  timechrx <- Sys.time()
  chr = obj$markers[1,2]
  
  setwd(outdir)
  
  file.prefix = paste(tx, pheno.col, sep = "_")
  
  plot.title = paste(tx, pheno.col, sep = " ")
  
  strains = sub("/", "_", hs.colors[,2])
  
  load(file = paste0(sanger.dir, "X.Rdata"))
  
  null.mod = glm(pheno[,pheno.col] ~ addcovar, family = binomial(logit))
  #null.mod = glm(trait ~ addcovar, family = poisson(link = "log"))
  null.ll = logLik(null.mod)
  pv = rep(0, nrow(sanger))
  
  glm.fxn = function(snp.rng, local.probs) {
    
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
    sdps.to.use = which(rowSums(cur.alleles) > 1.0)
    
    sex.col = which(colnames(addcovar) == "sex")
    if(length(sex.col) != 1) {
      stop("One of the columns of addcovar MUST be named 'sex'.")
    } # if(length(sex.col) != 1)
    
    # Run the model at each unique SDP.
    for(j in sdps.to.use) {
      
      
      full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = binomial(logit))
      #full.mod = glm(trait ~ addcovar + cur.alleles[j,], family = poisson(link = "log"))
      cur.ll[j] = logLik(full.mod)
      
    } # for(j)
    
    # This is the LRS.
    cur.ll = cur.ll - null.ll
    
    # Return the results.
    cur.ll[m]
    
  } # glm.fxn()
  
  # SNPs before the first marker.
  snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,1])
    
  } # if(length(snp.rng) > 0)
  
  # SNPs between Markers.
  for(i in 1:(nrow(obj$markers)-1)) {
    
    snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                      sanger.hdr$POS <= obj$markers[i+1,3])
    
    if(length(snp.rng) > 0) {
      
      # Take the mean of the haplotype probs at the surrounding markers.
      pv[snp.rng] = glm.fxn(snp.rng, (obj$probs[,,i] +
                                        obj$probs[,,i+1]) * 0.5)
      
    } # if(length(snp.rng) > 0)
    
  } # for(i)
  
  # SNPs after the last marker.
  snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])
    
  } # if(length(snp.rng) > 0)
  
  # Convert LRS to p-values using the chi-squared distribution.
  pv = pchisq(2 * pv, df = 1, lower.tail = FALSE)
  pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)
  
  maxPV = max(-log10(pv[,6]))
  print(paste0("Maximum LOD: ", maxPV), digits = 1)
  rm(maxPV)
  
  save(pv, file = paste0(file.prefix, "_chr", chr, ".Rdata"))
  
  png(paste0(file.prefix, "_chr", chr,".png"), width = 2600,
      height = 1200, res = 130)
  plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
  mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
  dev.off()
  
  # Return the positions and p-values.
  return(pv)
  rm(sanger, sanger.hdr)
  print(paste(round(difftime(Sys.time(), timechrx, units = 'hours'), digits = 2)))
  
} # GRSDbinom.xchr.fast()
################################################################################
synch.sample.IDs = function(pheno = NULL, probs = NULL, K = NULL, addcovar = NULL,
                            intcovar = NULL) {
  samples = NULL
  # pheno
  if(!is.null(pheno)) {
    samples = rownames(pheno)
  } # if(!is.null(pheno))
  # probs
  if(!is.null(probs)) {
    if(is.null(samples)) {
      samples = rownames(probs)
    } else {
      samples = intersect(samples, rownames(probs))
    } # else
  } # if(!is.null(probs))
  # K
  if(!is.null(K)) {
    if(is.null(samples)) {
      if(is.list(K)) {
        samples = rownames(K[[1]])
      } else {
        samples = rownames(K)
      } # else
    } else {
      if(is.list(K)) {
        samples = intersect(samples, rownames(K[[1]]))
      } else {
        samples = intersect(samples, rownames(K))
      } # else
    } # else
  } # if(!is.null(probs))
  # addcovar
  if(!is.null(addcovar)) {
    if(is.null(samples)) {
      samples = rownames(addcovar)
    } else {
      samples = intersect(samples, rownames(addcovar))
    } # else
  } # if(!is.null(probs))
  # intcovar
  if(!is.null(intcovar)) {
    if(is.null(samples)) {
      samples = rownames(intcovar)
    } else {
      samples = intersect(samples, rownames(intcovar))
    } # else
  } # if(!is.null(probs))
  if(length(samples) == 0) {
    warning(paste("There were no samples in common among the variables",
                  "passed in. Please make sure that there are rownames in",
                  "common between all variables."))
  } # if(length(samples) == 0)
  # Now subset all of the variables and return them in a list.
  retval = NULL
  if(!is.null(pheno)) {
    pheno = pheno[samples,,drop = FALSE]
    retval = list(pheno = pheno)
  } # if(!is.null(pheno))
  if(!is.null(probs)) {
    probs = probs[samples,,,drop = FALSE]
    if(is.null(retval)) {
      retval = list(probs = probs)
    } else {
      retval[[length(retval) + 1]] = probs
      names(retval)[length(retval)] = "probs"
    } # else
  } # if(!is.null(probs))
  if(!is.null(K)) {
    if(is.list(K)) {
      for(i in 1:length(K)) {
        K[[i]] = K[[i]][samples, samples]
      } # for(i)
    } else {
      K = K[samples, samples]
    } # else
    if(is.null(retval)) {
      retval = list(K = K)
    } else {
      retval[[length(retval) + 1]] = K
      names(retval)[length(retval)] = "K"
    } # else
  } # if(!is.null(K))
  if(!is.null(addcovar)) {
    addcovar = addcovar[samples,,drop = FALSE]
    if(is.null(retval)) {
      retval = list(addcovar = addcovar)
    } else {
      retval[[length(retval) + 1]] = addcovar
      names(retval)[length(retval)] = "addcovar"
    } # else
  } # if(!is.null(addcovar))
  if(!is.null(intcovar)) {
    intcovar = intcovar[samples,,drop = FALSE]
    if(is.null(retval)) {
      retval = list(intcovar = intcovar)
    } else {
      retval[[length(retval) + 1]] = intcovar
      names(retval)[length(retval)] = "intcovar"
    } # else
  } # if(!is.null(intcovar))
  return(retval)
} # synch.sample.IDs()
################################################################################
scanone = function(pheno, pheno.col = 1, probs = NULL, K = NULL, addcovar = NULL,
                   intcovar = NULL, snps = NULL, model = c("additive", "full")) {
  model = match.arg(model)
  if(is.null(rownames(pheno))) {
    stop("rownames(pheno) is null. The sample IDs must be in rownames(pheno).")
  } # if(is.null(rownames(pheno)))
  if(is.null(addcovar)) {
    stop(paste("You must map using \\'sex\\' as an additive covariate. Also,",
               "we require sex to map on the X chromosome. We even require sex if",
               "you are only mapping with one sex."))
  } # if(missing(addcovar))
  sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
  if(length(sex.col) == 0) {
    stop(paste("addcovar must contain a column called \\'sex\\'. Please add",
               "a sex column with sex coded as either F & M or 0 for females and 1 for males."))
  } # if(length(sex.col) == 0)
  if(is.null(rownames(addcovar))) {
    stop("rownames(addcovar) is null. The sample IDs must be in rownames(addcovar).")
  } # if(is.null(rownames(addcovar)))
  if(!is.null(intcovar)) {
    if(is.null(rownames(intcovar))) {
      stop("rownames(intcovar) is null. The sample IDs must be in rownames(intcovar).")
    } # if(is.null(rownames(intcovar)))
  } # if(!is.null(intcovar))
  snps[,2] = as.character(snps[,2])
  num.auto = get.num.auto(snps)
  # Convert phenotype names to phenotype column numbers.
  if(is.character(pheno.col)) {
    pheno.col = match(pheno.col, colnames(pheno))
  } # if(is.character(pheno.col))
  # Get the intersection of all of the sample IDs from pheno, probs, K, 
  # addcovar and intcovar.
  tmp = synch.sample.IDs(pheno = pheno, probs = probs, K = K, addcovar = addcovar,
                         intcovar = intcovar)
  pheno = tmp$pheno
  probs = tmp$probs
  K = tmp$K
  addcovar = tmp$addcovar
  intcovar = tmp$intcovar
  print(paste("Mapping with", nrow(pheno),"samples."))
  if(any(dim(probs) == 0)) {
    stop(paste("There are no matching samples in the data. Please",
               "verify that the sample IDs in rownames(pheno) match the sample",
               "IDs in rownames(probs), rownames(addcovar) and rownames(K)."))
  } # if(any(dim(probs) == 0))
  snps = snps[snps[,1] %in% dimnames(probs)[[3]],]
  probs = probs[,,match(snps[,1], dimnames(probs)[[3]])]
  print(paste("Mapping with", nrow(snps), "markers."))
  if(any(dim(probs) == 0)) {
    stop(paste("There are no matching markers in snps and probs. Please",
               "verify that the marker IDs in snps[,1] match the marker",
               "IDs in dimnames(probs)[[3]]."))
  } # if(any(dim(probs) == 0))
  if(sum(rownames(pheno) %in% rownames(addcovar)) == 0) {
    stop(paste("rownames(pheno) does not contain any sample IDs in",
               "common with rownames(addcovar). Please make sure that the",
               "rownames in pheno and addcovar match."))
  } # if(sum(rownames(pheno) %in% rownames(addcovar)) == 0)
  addcovar = as.matrix(addcovar)
  # Match sample IDs in interactive covariates.
  if(!is.null(intcovar)) {
    intcovar = as.matrix(intcovar)
    intcovar = intcovar[rownames(intcovar) %in% rownames(pheno),,drop = FALSE]
    intcovar = intcovar[match(rownames(pheno), rownames(intcovar)),,drop = FALSE]
    if(is.null(colnames(intcovar))) {
      colnames(intcovar) = paste("intcovar", 1:ncol(intcovar), sep = ".")
    } # if(is.null(colnames(intcovar)))
  } # if(!is.null(intcovar))
  if(!is.null(K)) {
    # LOCO method.
    if(is.list(K)) {
      for(c in 1:length(K)) {
        K[[c]] = K[[c]][rownames(K[[c]]) %in% rownames(pheno), 
                        colnames(K[[c]]) %in% rownames(pheno)]
        K[[c]] = K[[c]][match(rownames(pheno), rownames(K[[c]])), 
                        match(rownames(pheno), colnames(K[[c]]))]
      } # for(c)
    } else {
      K = K[rownames(K) %in% rownames(pheno), colnames(K) %in% rownames(pheno)]
      K = K[match(rownames(pheno), rownames(K)), match(rownames(pheno), colnames(K))]
    } # else
  } # if(!is.null(K))
  retval = NULL
  if(is.null(K)) {
    retval = scanone.noK(pheno, pheno.col, probs, addcovar, intcovar, snps, model)
  } else if(is.list(K)) {
    retval = scanone.LOCO(pheno, pheno.col, probs, K, addcovar, intcovar, snps, model)
  } else {
    retval = scanone.K(pheno, pheno.col, probs, K, addcovar, intcovar, snps, model)
  } # else
  if(length(retval) == 1) {
    retval = retval[[1]]
  } # if(length(retval) == 1)
  return(retval)
} # scanone()
################################################################################
scanone.noK = function(pheno, pheno.col, probs, addcovar, intcovar, snps, model) {
  num.auto = get.num.auto(snps)
  xchr = which(snps[,2] %in% "X")
  # We require sex to be in addcovar.
  if(length(xchr) > 0) {
    sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
    if(length(sex.col) == 0) {
      stop(paste("You must map using \\'sex\\' as an additive covariate. Also,",
                 "we require sex to map on the X chromosome. We even require sex if",
                 "you are only mapping with one sex."))
    } # if(length(sex.col) == 0)
    addcovar[,sex.col] = as.numeric(factor(addcovar[,sex.col])) - 1
  } # if(length(xchr) > 0)
  # Make the results list.
  retval = as.list(1:length(pheno.col))
  names(retval) = colnames(pheno)[pheno.col]
  index = 1
  for(i in pheno.col) {
    print(colnames(pheno)[i])
    p = pheno[,i]
    names(p) = rownames(pheno)
    keep = which(!is.na(p) & !is.nan(p) & !is.infinite(p))
    # Autosomes
    auto.qtl = NULL
    if(!is.na(num.auto)) {
      auto = which(snps[,2] %in% 1:num.auto)
      # With covariates.
      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0 & 
                                     rowSums(is.nan(addcovar)) == 0 &
                                     rowSums(is.infinite(addcovar)) == 0))
      if(is.null(intcovar)) {
        # Additive covariates only.
        auto.qtl = fast.qtlrel(pheno = p[keep], probs = probs[keep,,auto], 
                               addcovar = addcovar[keep,,drop = FALSE],
                               snps = snps[auto,])
      } else {
        auto.qtl = qtl.qtlrel(pheno = p[keep], probs = probs[keep,,auto],
                              addcovar = addcovar[keep,,drop = FALSE], 
                              intcovar = intcovar[keep,,drop = FALSE], 
                              snps = snps[auto,])
      } # else
      auto.qtl = list(lod = list(A = auto.qtl$lod), 
                      coef = list(A = auto.qtl$coef))
    } # if(!is.na(num.auto))
    # X chromosome.
    if(length(xchr) > 0) {
      # Get the sex from addcovar. We forced it to be 0 or 1 above.
      sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
      females = which(addcovar[,sex.col] == 0)
      males   = which(addcovar[,sex.col] == 1)
      mfprobs = NULL
      if(length(females) > 0 & length(males) > 0) {
        # Take the male and female probabilities and place them into 
        # one big array.
        if(model == "additive") {
          mfprobs = array(0, c(dim(probs)[1], 2 * dim(probs)[2], length(xchr)),
                          dimnames = list(dimnames(probs)[[1]], paste(rep(c("F", "M"), 
                                                                          each = dim(probs)[2]), dimnames(probs)[[2]], sep = "."), 
                                          dimnames(probs)[[3]][xchr]))
          for(j in females) {
            mfprobs[j,1:dim(probs)[2],] = probs[j,,xchr]
          } # for(j)
          for(j in males) {
            mfprobs[j,(dim(probs)[2] + 1):dim(mfprobs)[2],] = probs[j,,xchr]
          } # for(j)
        } else if(model == "full") {
          tmp = matrix(unlist(strsplit(dimnames(probs)[[2]], split = "")),
                       nrow = 2)
          homo = dimnames(probs)[[2]][which(tmp[1,] == tmp[2,])]
          mfprobs = array(0, c(dim(probs)[1], dim(probs)[2] + length(homo),
                               length(xchr)),
                          dimnames = list(dimnames(probs)[[1]], c(paste(rep("F", 
                                                                            each = dim(probs)[2]), dimnames(probs)[[2]], sep = "."),
                                                                  paste("M", homo, sep = ".")), dimnames(probs)[[3]][xchr]))
          for(j in females) {
            mfprobs[j,1:dim(probs)[2],] = probs[j,,xchr]
          } # for(j)
          for(j in males) {
            mfprobs[j,(dim(probs)[2] + 1):dim(mfprobs)[2],] = probs[j,homo,xchr]
          } # for(j)
        } # else if(model == "full")
        # If we have both males and females, then we need to remove one 
        # column from the males.
        mfprobs = mfprobs[,-grep("M.A", dimnames(mfprobs)[[2]]),]
      } else if(length(females) > 0) {
        mfprobs = probs[,,xchr]
        dimnames(mfprobs)[[2]] = paste("F", dimnames(mfprobs)[[2]], sep = ".")
      } else if(length(males) > 0) {
        mfprobs = probs[,,xchr]
        dimnames(mfprobs)[[2]] = paste("M", dimnames(mfprobs)[[2]], sep = ".")
      } # else if(length(males) > 0)
      # With covariates.
      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))
      if(is.null(intcovar)) {
        # Additive covariates only.
        x.qtl = fast.qtlrel(pheno = p[keep], probs = mfprobs[keep,,], 
                            addcovar = addcovar[keep,-1,drop = FALSE],
                            snps = snps[xchr,])
      } else {
        # Additive & interactive covariates.
        x.qtl = qtl.qtlrel(pheno = p[keep], probs = mfprobs[keep,,],
                           addcovar = addcovar[keep,,drop = FALSE], 
                           intcovar = intcovar[keep,,drop = FALSE], snps = snps[xchr,])
      } # else
      if(!is.null(auto.qtl)) {
        auto.qtl$lod  = list(A = auto.qtl$lod$A,  X = x.qtl$lod)
        auto.qtl$coef = list(A = auto.qtl$coef$A, X = x.qtl$coef)
      } else {
        auto.qtl = x.qtl
      } # else
    } # if(length(xchr) > 0)
    retval[[index]] = auto.qtl
    class(retval[[index]]) = c("doqtl", class(retval[[index]]))
    attr(retval[[index]], "model") = "additive"
    index = index + 1
  } # for(i)
  if(length(retval) == 1) {
    retval = retval[[1]]
  } # if(length(retval) == 1)
  return(retval)
} # scanone.noK()
################################################################################
scanone.K = function(pheno, pheno.col = 1, probs, K, addcovar, intcovar, snps, model) {
  num.auto = get.num.auto(snps)
  xchr = which(snps[,2] %in% "X")
  # We require sex to be in addcovar.
  sex = NULL
  if(length(xchr) > 0) {
    sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
    if(length(sex.col) == 0) {
      stop(paste("You must map using \\'sex\\' as an additive covariate. Also,",
                 "we require sex to map on the X chromosome. We even require sex if",
                 "you are only mapping with one sex."))
    } # if(length(sex.col) == 0)
    addcovar[,sex.col] = as.numeric(factor(addcovar[,sex.col])) - 1
  } # if(length(xchr) > 0)
  # Make the results list.
  retval = as.list(1:length(pheno.col))
  names(retval) = colnames(pheno)[pheno.col]
  index = 1
  for(i in pheno.col) {
    print(colnames(pheno)[i])
    p = pheno[,i]
    names(p) = rownames(pheno)
    keep = which(!is.na(p) & !is.nan(p) & !is.infinite(p))
    # Autosomes
    auto.qtl = NULL
    if(!is.na(num.auto)) {
      auto = which(snps[,2] %in% 1:num.auto)
      # With covariates.
      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0 & 
                                     rowSums(is.nan(addcovar)) == 0 &
                                     rowSums(is.infinite(addcovar)) == 0))
      if(is.null(intcovar)) {
        # Additive covariates only.
        auto.qtl = fast.qtlrel(pheno = p[keep], probs = probs[keep,,auto], 
                               K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                               snps = snps[auto,])
      } else {
        auto.qtl = qtl.qtlrel(pheno = p[keep], probs = probs[keep,,auto],
                              K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE], 
                              intcovar = intcovar[keep,,drop = FALSE], snps = snps[auto,])
      } # else
      auto.qtl = list(lod = list(A = auto.qtl$lod), 
                      coef = list(A = auto.qtl$coef))
    } # if(!is.na(num.auto))
    # X chromosome.
    if(length(xchr) > 0) {
      # Get the sex from addcovar. We forced it to be 0 or 1 above.
      sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
      females = which(addcovar[,sex.col] == 0)
      males   = which(addcovar[,sex.col] == 1)
      mfprobs = NULL
      if(length(females) > 0 & length(males) > 0) {
        # Take the male and female probabilities and place them into 
        # one big array.
        if(model == "additive") {
          mfprobs = array(0, c(dim(probs)[1], 2 * dim(probs)[2], length(xchr)),
                          dimnames = list(dimnames(probs)[[1]], paste(rep(c("F", "M"), 
                                                                          each = dim(probs)[2]), dimnames(probs)[[2]], sep = "."), 
                                          dimnames(probs)[[3]][xchr]))
          for(j in females) {
            mfprobs[j,1:dim(probs)[2],] = probs[j,,xchr]
          } # for(j)
          for(j in males) {
            mfprobs[j,(dim(probs)[2] + 1):dim(mfprobs)[2],] = probs[j,,xchr]
          } # for(j)
        } else if(model == "full") {
          tmp = matrix(unlist(strsplit(dimnames(probs)[[2]], split = "")),
                       nrow = 2)
          homo = dimnames(probs)[[2]][which(tmp[1,] == tmp[2,])]
          mfprobs = array(0, c(dim(probs)[1], dim(probs)[2] + length(homo),
                               length(xchr)),
                          dimnames = list(dimnames(probs)[[1]], c(paste(rep("F", 
                                                                            each = dim(probs)[2]), dimnames(probs)[[2]], sep = "."),
                                                                  paste("M", homo, sep = ".")), dimnames(probs)[[3]][xchr]))
          for(j in females) {
            mfprobs[j,1:dim(probs)[2],] = probs[j,,xchr]
          } # for(j)
          for(j in males) {
            mfprobs[j,(dim(probs)[2] + 1):dim(mfprobs)[2],] = probs[j,homo,xchr]
          } # for(j)
        } # else if(model == "full")
        # If we have both males and females, then we need to remove one 
        # column from the males.
        mfprobs = mfprobs[,-grep("M.A", dimnames(mfprobs)[[2]]),]
      } else if(length(females) > 0) {
        mfprobs = probs[,,xchr]
        dimnames(mfprobs)[[2]] = paste("F", dimnames(mfprobs)[[2]], sep = ".")
      } else if(length(males) > 0) {
        mfprobs = probs[,,xchr]
        dimnames(mfprobs)[[2]] = paste("M", dimnames(mfprobs)[[2]], sep = ".")
      } # else if(length(males) > 0)
      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0))
      if(is.null(intcovar)) {
        # Additive covariates only.
        x.qtl = fast.qtlrel(pheno = p[keep], probs = mfprobs[keep,,], 
                            K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                            snps = snps[xchr,])
      } else {
        # Additive & interactive covariates.
        x.qtl = qtl.qtlrel(pheno = p[keep], probs = mfprobs[keep,,],
                           K = K[keep,keep], addcovar = addcovar[keep,,drop = FALSE], 
                           intcovar = intcovar[keep,,drop = FALSE], snps = snps[xchr,])
      } # else
      if(!is.null(auto.qtl)) {
        auto.qtl$lod  = list(A = auto.qtl$lod$A,  X = x.qtl$lod)
        auto.qtl$coef = list(A = auto.qtl$coef$A, X = x.qtl$coef)
      } else {
        auto.qtl = x.qtl
      } # else
    } # if(length(xchr) > 0)
    retval[[index]] = auto.qtl
    class(retval[[index]]) = c("doqtl", class(retval[[index]]))
    attr(retval[[index]], "model") = "additive"
    index = index + 1
  } # for(i)
  if(length(retval) == 1) {
    retval = retval[[1]]
  } # if(length(retval) == 1)
  return(retval)
} # scanone.K()
################################################################################
scanone.LOCO = function(pheno, pheno.col = 1, probs, K, addcovar, intcovar, snps, model) {
  num.auto = get.num.auto(snps)
  xchr = which(snps[,2] %in% "X")
  # We require sex to be in addcovar.
  sex = NULL
  if(length(xchr) > 0) {
    sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
    if(length(sex.col) == 0) {
      stop(paste("You must map using \\'sex\\' as an additive covariate. Also,",
                 "we require sex to map on the X chromosome. We even require sex if",
                 "you are only mapping with one sex."))
    } # if(length(sex.col) == 0)
    addcovar[,sex.col] = as.numeric(factor(addcovar[,sex.col])) - 1
  } # if(length(xchr) > 0)
  # Make the results list.
  retval = as.list(1:length(pheno.col))
  names(retval) = colnames(pheno)[pheno.col]
  index = 1
  for(i in pheno.col) {
    print(colnames(pheno)[i])
    p = pheno[,i]
    names(p) = rownames(pheno)
    keep = which(!is.na(p) & !is.nan(p) & !is.infinite(p))
    # Autosomes
    auto.qtl = NULL
    if(!is.na(num.auto)) {
      keep = intersect(keep, which(rowSums(is.na(addcovar)) == 0 & 
                                     rowSums(is.nan(addcovar)) == 0 &
                                     rowSums(is.infinite(addcovar)) == 0))
      if(is.null(intcovar)) {
        # Additive covariates only.
        snprng = which(snps[,2] == 1)
        auto.qtl = fast.qtlrel(pheno = p[keep], probs = probs[keep,,snprng], 
                               K = K[[1]][keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                               snps = snps[snprng,])
        for(c in 2:num.auto) {
          snprng = which(snps[,2] == c)
          tmp = fast.qtlrel(pheno = p[keep, drop = FALSE], probs = probs[keep,,snprng], 
                            K = K[[c]][keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                            snps = snps[snprng,])
          auto.qtl$lod  = rbind(auto.qtl$lod,  tmp$lod)
          auto.qtl$coef = rbind(auto.qtl$coef, tmp$coef)
        } # for(c)
      } else {
        # Additive and interactive covariates.
        snprng = which(snps[,2] == 1)
        auto.qtl = qtl.qtlrel(pheno = p[keep, drop = FALSE], 
                              probs = probs[keep,,snprng], K = K[[1]][keep,keep], 
                              addcovar = addcovar[keep,,drop = FALSE], 
                              intcovar = intcovar[keep,,drop = FALSE], snps = snps[snprng,])
        for(c in 2:num.auto) {
          snprng = which(snps[,2] == c)
          tmp = qtl.qtlrel(pheno = p[keep, drop = FALSE],
                           probs = probs[keep,,snprng], K = K[[c]][keep,keep], 
                           addcovar = addcovar[keep,,drop = FALSE], 
                           intcovar = intcovar[keep,,drop = FALSE], snps = snps[snprng,])
          auto.qtl$lod  = rbind(auto.qtl$lod,  tmp$lod)
          auto.qtl$coef = rbind(auto.qtl$coef, tmp$coef)
        } # for(c)
      } # else
      auto.qtl = list(lod  = list(A = auto.qtl$lod), 
                      coef = list(A = auto.qtl$coef))
    } # if(!is.na(num.auto))
    # X chromosome.
    if(length(xchr) > 0) {
      # Get the sex from addcovar. We forced it to be 0 or 1 above.
      sex.col = grep("^sex$", colnames(addcovar), ignore.case = TRUE)
      females = which(addcovar[,sex.col] == 0)
      males   = which(addcovar[,sex.col] == 1)
      mfprobs = NULL
      if(length(females) > 0 & length(males) > 0) {
        # Take the male and female probabilities and place them into 
        # one big array.
        if(model == "additive") {
          mfprobs = array(0, c(nrow(probs), 2 * ncol(probs), length(xchr)),
                          dimnames = list(dimnames(probs)[[1]], paste(rep(c("F", "M"), 
                                                                          each = ncol(probs)), dimnames(probs)[[2]], sep = "."), 
                                          dimnames(probs)[[3]][xchr]))
          for(j in females) {
            mfprobs[j,1:dim(probs)[2],] = probs[j,,xchr]
          } # for(j)
          for(j in males) {
            mfprobs[j,(dim(probs)[2] + 1):dim(mfprobs)[2],] = probs[j,,xchr]
          } # for(j)
        } else if(model == "full") {
          tmp = matrix(unlist(strsplit(dimnames(probs)[[2]], split = "")),
                       nrow = 2)
          homo = dimnames(probs)[[2]][which(tmp[1,] == tmp[2,])]
          mfprobs = array(0, c(dim(probs)[1], dim(probs)[2] + length(homo),
                               length(xchr)),
                          dimnames = list(dimnames(probs)[[1]], c(paste(rep("F", 
                                                                            each = dim(probs)[2]), dimnames(probs)[[2]], sep = "."),
                                                                  paste("M", homo, sep = ".")), dimnames(probs)[[3]][xchr]))
          for(j in females) {
            mfprobs[j,1:dim(probs)[2],] = probs[j,,xchr]
          } # for(j)
          for(j in males) {
            mfprobs[j,(dim(probs)[2] + 1):dim(mfprobs)[2],] = probs[j,homo,xchr]
          } # for(j)
        } # else if(model == "full")
        # If we have both males and females, then we need to remove one 
        # column from the males.
        mfprobs = mfprobs[,-grep("M.A", dimnames(mfprobs)[[2]]),]
      } else if(length(females) > 0) {
        mfprobs = probs[,,xchr]
        dimnames(mfprobs)[[2]] = paste("F", dimnames(mfprobs)[[2]], sep = ".")
      } else if(length(males) > 0) {
        mfprobs = probs[,,xchr]
        dimnames(mfprobs)[[2]] = paste("M", dimnames(mfprobs)[[2]], sep = ".")
      } # else if(length(males) > 0)
      x.qtl = NULL
      if(is.null(intcovar)) {
        # Additive covariates only.
        x.qtl = qtl.qtlrel(pheno = p[keep], probs = mfprobs[keep,,], 
                           K = K[["X"]][keep,keep], addcovar = addcovar[keep,,drop = FALSE],
                           snps = snps[xchr,])
      } else {
        # Additive & interactive covariates.
        x.qtl = qtl.qtlrel(pheno = p[keep], probs = mfprobs[keep,,],
                           K = K[["X"]][keep,keep], addcovar = addcovar[keep,,drop = FALSE], 
                           intcovar = intcovar[keep,,drop = FALSE], snps = snps[xchr,])
      } # else
      if(!is.null(auto.qtl)) {
        auto.qtl$lod  = list(A = auto.qtl$lod$A,  X = x.qtl$lod)
        auto.qtl$coef = list(A = auto.qtl$coef$A, X = x.qtl$coef)
      } else {
        auto.qtl = x.qtl
      } # else
    } # if(length(xchr) > 0)
    retval[[index]] = auto.qtl
    class(retval[[index]]) = c("doqtl", class(retval[[index]]))
    attr(retval[[index]], "model") = "additive"
    index = index + 1
  } # for(i)
  if(length(retval) == 1) {
    retval = retval[[1]]
  } # if(length(retval) == 1)
  return(retval)
} # scanone.LOCO()
################################################################################
get.num.auto = function(snps) {
  # Turn off the warnings for a moment while we get the number of autosomes.
  old.warn = options("warn")$warn
  options(warn = -1)
  num.auto = max(as.numeric(snps[,2]), na.rm = TRUE)
  options(warn = old.warn)
  return(num.auto)
  
} # get.num.auto()
################################################################################
matrixeqtl.snps = function(pheno, geno, K, addcovar) {
  pheno = as.matrix(t(pheno))
  geno  = as.matrix(t(geno))
  
  # Create an error covariance matrix.
  if(!missing(K)) {
    eig = eigen(K, symmetric = TRUE)
    if(any(eig$values <= 0)) {
      stop("The covariance matrix is not positive definite")
    } # if(any(eig$values <= 0))
    correctionMatrix = eig$vectors %*% diag(1.0 / sqrt(eig$values)) %*% 
      t(eig$vectors)
    rm(eig)
  } else {
    correctionMatrix = numeric()
  } # else
  
  # Add an intercept and rotate the covariates.
  cvrt = matrix(1, nrow = 1, ncol = ncol(pheno))
  if(!missing(addcovar)) {
    cvrt = rbind(cvrt, t(addcovar))
  } # if(!missing(addcovar))
  
  if(length(correctionMatrix) > 0) {
    cvrt = cvrt %*% correctionMatrix
  } # if(length(correctionMatrix) > 0)
  q = qr(t(cvrt))
  cvrt = t(qr.Q(q))
  
  # Rotate and center the genes.
  if(length(correctionMatrix) > 0) {
    pheno = pheno %*% correctionMatrix
  } # if(length(correctionMatrix) > 0)
  pheno = pheno - tcrossprod(pheno, cvrt) %*% cvrt
  div = sqrt(rowSums(pheno^2))
  div[div == 0] = 1
  pheno = pheno / div
  rm(div)
  # Rotate and center the SNPs. 
  if(length(correctionMatrix) > 0) {
    geno = geno %*% correctionMatrix
  } # if(length(correctionMatrix) > 0)
  geno = geno - tcrossprod(geno, cvrt) %*% cvrt
  div = sqrt(rowSums(geno^2))
  drop = div < 5 * sqrt(ncol(geno) * .Machine$double.eps)
  div[drop] = 1
  geno = geno / div
  geno[drop,] = 0
  # Note: we return the R^2.
  return(tcrossprod(geno, pheno)^2)
} # matrixeqtl.snps()
################################################################################
bayesint = function(qtl, chr, prob = 0.95, expandtomarkers = TRUE) {
  if(missing(qtl)) {
    stop("bayesint: The qtl cannot be null. Please supply a QTL object.")
  } # if(is.null(qtl))
  
  if(missing(chr)) {
    stop(paste("bayesint: The chromosome cannot be null."))
  } else if(!chr %in% c(1:19, "X")) {
    stop(paste("bayesint: The chromosome must be 1 to 19 or X."))
  } # else if
  
  if(prob < 0 | prob > 1) {
    stop(paste("bayesint: The probability must between 0 and 1."))
  } # if(prob < 0 | prob > 1)
  
  old.warn = options("warn")$warn
  options(warn = -1)
  if(!is.na(as.numeric(chr))) {
    qtl = qtl$lod$A
  } else {
    qtl = qtl$lod$X
  } # else
  options(warn = old.warn)
  
  qtl[,1] = as.character(qtl[,1])
  qtl[,2] = as.character(qtl[,2])
  qtl[,3] = as.numeric(qtl[,3])
  qtl[,7] = as.numeric(qtl[,7])
  qtl = qtl[qtl[,2] == chr,]
  pos = qtl[,3]
  
  if(any(is.na(pos))) {
    remove = which(is.na(pos))
    qtl = qtl[-remove,]
    pos = pos[-remove]
  } # if(any(is.na(pos)))
  
  # Make a set of 10,000 intervals so that we can integrate numerically.
  breaks = approx(x = pos, y = 10^qtl[,7], xout = seq(pos[1], pos[length(pos)],
                                                      length.out = 1e5))
  widths  = diff(breaks$x)
  heights = breaks$y[-1] + breaks$y[-length(breaks$y)]
  trapezoids = 0.5 * heights * widths
  # Normalize the area to 1.
  trapezoids = trapezoids / sum(trapezoids)
  # This code copied from the R/qtl bayesint() function by Karl Broman.
  ord = order(breaks$y[-length(breaks$y)], decreasing = TRUE)
  wh  = min(which(cumsum(trapezoids[ord]) >= prob))
  int = range(ord[1:wh])
  # Find the left & right SNP.
  left.snp  = c(NA, qtl[1,2], breaks$x[int][1], 
                approx(qtl[,3], qtl[,4], breaks$x[int][1])$y,
                approx(qtl[,3], qtl[,5], breaks$x[int][1])$y,
                approx(qtl[,3], qtl[,6], breaks$x[int][1])$y,
                approx(qtl[,3], qtl[,7], breaks$x[int][1])$y)
  max.snp   = qtl[which.max(qtl[,7]),]
  right.snp = c(NA, qtl[1,2], breaks$x[int][2], 
                approx(qtl[,3], qtl[,4], breaks$x[int][2])$y,
                approx(qtl[,3], qtl[,5], breaks$x[int][2])$y,
                approx(qtl[,3], qtl[,6], breaks$x[int][2])$y,
                approx(qtl[,3], qtl[,7], breaks$x[int][2])$y)
  if(expandtomarkers) {
    # Find the left & right SNP.
    left.snp  = qtl[max(which(breaks$x[int][1] >= qtl[,3])),]
    max.snp   = qtl[which.max(qtl[,7]),]
    right.snp = qtl[min(which(breaks$x[int][2] <= qtl[,3])),] 
  } # if(expandtomarkers)
  retval = rbind(left.snp, max.snp, right.snp)
  retval[,3] = round(as.numeric(retval[,3]), digits = 6)
  retval[,4] = round(as.numeric(retval[,4]), digits = 6)
  retval[,5] = round(as.numeric(retval[,5]), digits = 6)
  retval[,6] = round(as.numeric(retval[,6]), digits = 6)
  retval$lod = as.numeric(retval[,7])
  return(retval)
} # bayesint()
################################################################################
GRSD.poisson = function(pheno, pheno.col, probs, K, addcovar, markers, snp.file,
                        outdir = "~/Desktop/", tx = c("Gamma", "HZE", "Unirradiated"), 
                        sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/"){
  begin <- Sys.time()
  # COVARIATES #
  
  samples = intersect(rownames(pheno), rownames(probs))
  samples = intersect(samples, rownames(addcovar))
  samples = intersect(samples, rownames(K[[1]]))
  stopifnot(length(samples) > 0)
  print(paste("A total of", length(samples), tx, "samples are complete."))
  
  pheno = pheno[samples,,drop = FALSE]
  addcovar = addcovar[samples,,drop = FALSE]
  probs = probs[samples,,,drop = FALSE]
  
  
  
  # DEFINE TRAIT #
  
  file.prefix = paste(tx, pheno.col, "poisson", sep = "_")
  
  plot.title = paste(tx, pheno.col, sep = " ")
  print(plot.title)
  
  trait = pheno[,pheno.col]
  print(table(trait))
  
  
  # LOGISTIC REGRESSION MODEL #
  
  for(i in 1:length(K)) {
    K[[i]] = K[[i]][samples, samples]
  } # for(i)
  
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
  
  sanger.dir = sanger.dir
  
  # MAPPING ANALYSES #
  
  result = vector("list", length(data))
  names(result) = names(data)
  print(paste("Mapping with", length(samples), tx, "samples."))
  
  for(i in 1:19) {
    print(i)
    begin1 <- Sys.time()
    result[[i]] = GRSDpoisson(data[[i]], pheno, pheno.col, addcovar, tx, sanger.dir)
    print(paste(round(difftime(Sys.time(), begin1, units = 'mins'), digits = 2),
                "minutes..."))
  } #for(i)
  
  print("X")
  result[["X"]] = GRSDpoisson.xchr(data[["X"]], pheno, pheno.col, addcovar, tx, sanger.dir)
  
  print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
              "hours elapsed during mapping."))
  
  
  # PLOTTING
  plotter <- Sys.time()
  
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
  print("Plotting...")
  for(i in 1:length(files)) {
    
    load(files[i])
    data[[i]] = pv
    data[[i]][,6] = -log10(data[[i]][,6])
    
  } # for(i)
  
  num.snps = sapply(data, nrow)
  chrs = c(1:19, "X")
  
  xlim = c(0, sum(num.snps))
  ylim = c(0, max(sapply(data, function(z) { max(z[,6]) })))
  
  # PLOT ALL CHROMOSOMES #
  setwd(outdir)
  chrlen = get.chr.lengths()[1:20]
  chrsum = cumsum(chrlen)
  chrmid = c(1, chrsum[-length(chrsum)]) + chrlen * 0.5
  names(chrmid) = names(chrlen)
  
  png(paste0(file.prefix, "_QTL.png"), width = 2600, height = 1200, res = 200)
  plot(-1, -1, col = 0, xlim = c(0, max(chrsum)), ylim = ylim, xlab = "",
       ylab = "-log10(p-value)", las = 1, main = plot.title, xaxt = "n")
  for(i in 1:length(data)) {
    print(paste("Plotting chromosome", i))
    pos = data[[i]][,3] * 1e-6 + c(0, chrsum)[i]
    points(pos, data[[i]][,6], col = c("black", "grey50")[i %% 2 + 1],
           pch = 20)
  } # for(i)
  mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.5)
  dev.off()
  
  # Convert to GRangesList for storage
  chrs = c(1:19, "X")
  qtl = GRangesList(GRanges("list", length(result)))
  
  for(i in 1:length(chrs)) {
    print(i)
    qtl[[i]] <- GRanges(seqnames = Rle(result[[i]]$CHR),
                        ranges = IRanges(start = result[[i]]$POS, width = 1),
                        p.value = result[[i]]$pv)
  } # for(i)
  
  save(qtl, file.prefix, file = paste0(file.prefix, "_QTL.Rdata"))
  
  print(paste(round(difftime(Sys.time(), plotter, units = 'hours'), digits = 2),
              "hours elapsed during plotting."))
  
  
}
################################################################################
GRSDpoisson = function(obj, pheno, pheno.col, addcovar, tx, sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/") {
  chr = obj$markers[1,2]
  
  setwd(outdir)
  
  file.prefix = paste(tx, pheno.col, "poisson", sep = "_")
  
  plot.title = paste(tx, pheno.col, sep = " ")
  
  strains = sub("/", "_", hs.colors[,2])
  
  load(file = paste0(sanger.dir, chr, ".Rdata"))
  
  #null.mod = glm(pheno[,pheno.col] ~ addcovar, family = binomial(logit))
  null.mod = glm(pheno[,pheno.col] ~ addcovar, family = poisson(link = "log"))
  null.ll = logLik(null.mod)
  pv = rep(0, nrow(sanger))
  
  glm.fxn = function(snp.rng, local.probs) {
    
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
    sdps.to.use = which(rowSums(cur.alleles) > 1.0)
    
    # Run the model at each unique SDP.
    for(j in sdps.to.use) {
      
      
      #full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = binomial(logit))
      full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = poisson(link = "log"))
      cur.ll[j] = logLik(full.mod)
      
    } # for(j)
    
    # This is the LRS.
    cur.ll = cur.ll - null.ll
    
    # Return the results.
    cur.ll[m]
    
  } # glm.fxn()
  
  # SNPs before the first marker.
  snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,1])
    
  } # if(length(snp.rng) > 0)
  
  # SNPs between Markers.
  for(i in 1:(nrow(obj$markers)-1)) {
    
    snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                      sanger.hdr$POS <= obj$markers[i+1,3])
    
    if(length(snp.rng) > 0) {
      
      # Take the mean of the haplotype probs at the surrounding markers.
      pv[snp.rng] = glm.fxn(snp.rng, (obj$probs[,,i] +
                                        obj$probs[,,i+1]) * 0.5)
      
    } # if(length(snp.rng) > 0)
    
  } # for(i)
  
  # SNPs after the last marker.
  snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])
    
  } # if(length(snp.rng) > 0)
  
  # Convert LRS to p-values using the chi-squared distribution.
  pv = pchisq(1 * pv, df = 4, lower.tail = FALSE)
  pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)
  
  maxPV = max(-log10(pv[,6]))
  print(paste0("Maximum LOD: ", maxPV), digits = 1)
  rm(maxPV)
  
  save(pv, file = paste0(file.prefix, "_chr", chr, ".Rdata"))
  
  png(paste0(file.prefix, "_chr", chr,".png"), width = 2000,
      height = 1600, res = 200)
  plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
  mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
  dev.off()
  
  # Return the positions and p-values.
  return(pv)
  
} # GRSDpoisson()
################################################################################
GRSDpoisson.xchr = function(obj, pheno, pheno.col, addcovar, tx, sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/") {
  chr = obj$markers[1,2]
  
  setwd(outdir)
  
  file.prefix = paste(tx, pheno.col, "poisson", sep = "_")
  
  plot.title = paste(tx, pheno.col, sep = " ")
  
  strains = sub("/", "_", hs.colors[,2])
  
  load(file = paste0(sanger.dir, "X.Rdata"))
  
  #null.mod = glm(pheno[,pheno.col] ~ addcovar, family = binomial(logit))
  null.mod = glm(pheno[,pheno.col] ~ addcovar, family = poisson(link = "log"))
  null.ll = logLik(null.mod)
  pv = rep(0, nrow(sanger))
  
  glm.fxn = function(snp.rng, local.probs) {
    
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
    sdps.to.use = which(rowSums(cur.alleles) > 1.0)
    
    sex.col = which(colnames(addcovar) == "sex")
    if(length(sex.col) != 1) {
      stop("One of the columns of addcovar MUST be named 'sex'.")
    } # if(length(sex.col) != 1)
    
    # Run the model at each unique SDP.
    for(j in sdps.to.use) {
      
      
      #full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = binomial(logit))
      full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = poisson(link = "log"))
      cur.ll[j] = logLik(full.mod)
      
    } # for(j)
    
    # This is the LRS.
    cur.ll = cur.ll - null.ll
    
    # Return the results.
    cur.ll[m]
    
  } # glm.fxn()
  
  # SNPs before the first marker.
  snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,1])
    
  } # if(length(snp.rng) > 0)
  
  # SNPs between Markers.
  for(i in 1:(nrow(obj$markers)-1)) {
    
    snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                      sanger.hdr$POS <= obj$markers[i+1,3])
    
    if(length(snp.rng) > 0) {
      
      # Take the mean of the haplotype probs at the surrounding markers.
      pv[snp.rng] = glm.fxn(snp.rng, (obj$probs[,,i] +
                                        obj$probs[,,i+1]) * 0.5)
      
    } # if(length(snp.rng) > 0)
    
  } # for(i)
  
  # SNPs after the last marker.
  snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])
    
  } # if(length(snp.rng) > 0)
  
  # Convert LRS to p-values using the chi-squared distribution.
  pv = pchisq(1 * pv, df = 4, lower.tail = FALSE)
  pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)
  
  maxPV = max(-log10(pv[,6]))
  print(paste0("Maximum LOD: ", maxPV), digits = 1)
  rm(maxPV)
  
  save(pv, file = paste0(file.prefix, "_chr", chr, ".Rdata"))
  
  png(paste0(file.prefix, "_chr", chr,".png"), width = 2600,
      height = 1600, res = 200)
  plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
  mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
  dev.off()
  
  # Return the positions and p-values.
  return(pv)
  
} # GRSDpoisson.xchr()
################################################################################
GRSD.gauss = function(pheno, pheno.col, probs, K, addcovar, markers, snp.file,
                      outdir = "~/Desktop/", tx = c("Gamma", "HZE", "Unirradiated", "All"),
                      sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/"){
  begin <- Sys.time()
  begin
  # COVARIATES #
  
  samples = intersect(rownames(pheno), rownames(probs))
  samples = intersect(samples, rownames(addcovar))
  samples = intersect(samples, rownames(K[[1]]))
  stopifnot(length(samples) > 0)
  print(paste("A total of", length(samples), tx, "samples are complete."))
  
  pheno = pheno[samples,,drop = FALSE]
  addcovar = addcovar[samples,,drop = FALSE]
  probs = probs[samples,,,drop = FALSE]
  
  
  
  # DEFINE TRAIT #
  
  file.prefix = paste(tx, pheno.col, sep = "_")
  
  plot.title = paste(tx, pheno.col, sep = " ")
  print(plot.title)
  
  trait = pheno[,pheno.col]
  #print(table(trait))
  # print(paste(round(100*(sum(trait) / length(samples)), digits = 1),
  #             "% display the", pheno.col, "phenotype in the", tx, "group."))
  # 
  # LOGISTIC REGRESSION MODEL #
  
  for(i in 1:length(K)) {
    K[[i]] = K[[i]][samples, samples]
  } # for(i)
  
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
  
  setwd(outdir)
  
  # MAPPING ANALYSES #
  
  result = vector("list", length(data))
  names(result) = names(data)
  print(paste("Mapping with", length(samples), tx, "samples..."))
  sanger.dir = sanger.dir
  
  for(i in 1:19) {
    print(paste("CHROMOSOME", i))
    timechr <- Sys.time()
    result[[i]] = GRSDgauss.fast(data[[i]], pheno, pheno.col, addcovar, tx, sanger.dir)
    print(paste(round(difftime(Sys.time(), timechr, units = 'mins'), digits = 2),
                "minutes..."))
  } #for(i)
  
  print("X CHROMOSOME")
  result[["X"]] = GRSDgauss.xchr.fast(data[["X"]], pheno, pheno.col, addcovar, tx, sanger.dir)
  
  print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
              "hours elapsed during mapping."))
  
  # Convert to GRangesList for storage
  chrs = c(1:19, "X")
  qtl = GRangesList(GRanges("list", length(result)))
  
  for(i in 1:length(chrs)) {
    print(i)
    qtl[[i]] <- GRanges(seqnames = Rle(result[[i]]$CHR),
                        ranges = IRanges(start = result[[i]]$POS, width = 1),
                        p.value = result[[i]]$pv)
  } # for(i)
  
  
  
  # PLOTTING
  plotter <- Sys.time()
  
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
  print("Plotting...")
  for(i in 1:length(files)) {
    
    load(files[i])
    data[[i]] = pv
    data[[i]][,6] = -log10(data[[i]][,6])
    
  } # for(i)
  
  num.snps = sapply(data, nrow)
  chrs = c(1:19, "X")
  
  xlim = c(0, sum(num.snps))
  ylim = c(0, max(sapply(data, function(z) { max(z[,6]) })))
  
  # PLOT ALL CHROMOSOMES #
  setwd(outdir)
  chrlen = get.chr.lengths()[1:20]
  chrsum = cumsum(chrlen)
  chrmid = c(1, chrsum[-length(chrsum)]) + chrlen * 0.5
  names(chrmid) = names(chrlen)
  
  png(paste0(file.prefix, "_QTL.png"), width = 2600, height = 1200, res = 200)
  plot(-1, -1, col = 0, xlim = c(0, max(chrsum)), ylim = ylim, xlab = "",
       ylab = "-log10(p-value)", las = 1, main = plot.title, xaxt = "n")
  for(i in 1:length(data)) {
    print(paste("Plotting chromosome", i))
    pos = data[[i]][,3] * 1e-6 + c(0, chrsum)[i]
    points(pos, data[[i]][,6], col = c("black", "grey50")[i %% 2 + 1],
           pch = 20)
  } # for(i)
  mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.5)
  dev.off()
  
  # Convert to GRangesList for storage
  chrs = c(1:19, "X")
  qtl = GRangesList(GRanges("list", length(result)))
  
  for(i in 1:length(chrs)) {
    print(i)
    qtl[[i]] <- GRanges(seqnames = Rle(result[[i]]$CHR),
                        ranges = IRanges(start = result[[i]]$POS, width = 1),
                        p.value = result[[i]]$pv)
  } # for(i)
  
  save(qtl, file.prefix, file = paste0(file.prefix, "_QTL.Rdata"))
  
  print(paste(round(difftime(Sys.time(), plotter, units = 'hours'), digits = 2),
              "hours elapsed during plotting."))
  
}
################################################################################
GRSDgauss.fast = function(obj, pheno, pheno.col, addcovar, tx, sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/") {
  chr = obj$markers[1,2]
  
  setwd(outdir)
  
  file.prefix = paste(tx, pheno.col, sep = "_")
  
  plot.title = paste(tx, pheno.col, sep = " ")
  
  strains = sub("/", "_", hs.colors[,2])
  
  load(file = paste0(sanger.dir, chr, ".Rdata"))
  
  null.mod = glm(pheno[,pheno.col] ~ addcovar, family = gaussian)
  #null.mod = glm(trait ~ addcovar, family = poisson(link = "log"))
  null.ll = logLik(null.mod)
  pv = rep(0, nrow(sanger))
  
  glm.fxn = function(snp.rng, local.probs) {
    
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
    sdps.to.use = which(rowSums(cur.alleles) > 1.0)
    
    # Run the model at each unique SDP.
    for(j in sdps.to.use) {
      
      
      # library(regress)
      # mod = regress(obj$pheno[,obj$pheno.col] ~ obj$addcovar, ~obj$K, pos = c(TRUE, TRUE))
      # obj$K = mod$sigma[1] * obj$K + mod$sigma[2] * diag(nrow(obj$K))
      # rm(mod)
      
      
      
      full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = gaussian)
      #full.mod = glm(trait ~ addcovar + cur.alleles[j,], family = poisson(link = "log"))
      cur.ll[j] = logLik(full.mod)
      
    } # for(j)
    
    # This is the LRS.
    cur.ll = cur.ll - null.ll
    
    # Return the results.
    cur.ll[m]
    
  } # glm.fxn()
  
  # SNPs before the first marker.
  snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,1])
    
  } # if(length(snp.rng) > 0)
  
  # SNPs between Markers.
  for(i in 1:(nrow(obj$markers)-1)) {
    
    snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                      sanger.hdr$POS <= obj$markers[i+1,3])
    
    if(length(snp.rng) > 0) {
      
      # Take the mean of the haplotype probs at the surrounding markers.
      pv[snp.rng] = glm.fxn(snp.rng, (obj$probs[,,i] +
                                        obj$probs[,,i+1]) * 0.5)
      
    } # if(length(snp.rng) > 0)
    
  } # for(i)
  
  # SNPs after the last marker.
  snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])
    
  } # if(length(snp.rng) > 0)
  
  # Convert LRS to p-values using the chi-squared distribution.
  pv = pchisq(2 * pv, df = 1, lower.tail = FALSE)
  pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)
  
  max.LOD = max(-log10(pv[,6]))
  #max.LOD.position = pv[,2][which(max(-log10(pv[,6])))]
  
  print(paste0("Maximum LOD is ", max.LOD), digits = 1)
  rm(max.LOD, max.LOD.position)
  
  # save(pv, file = paste0(file.prefix, "_chr", chr, ".Rdata"))
  # 
  # png(paste0(file.prefix, "_chr", chr,".png"), width = 2600,
  #     height = 1200, res = 130)
  # plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
  # mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
  # dev.off()
  
  # Return the positions and p-values.
  return(pv)
  
  
} # GRSDgauss.fast()
################################################################################
GRSDgauss.xchr.fast = function(obj, pheno, pheno.col, addcovar, tx, sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/") {
  timechrx <- Sys.time()
  chr = obj$markers[1,2]
  
  setwd(outdir)
  
  file.prefix = paste(tx, pheno.col, sep = "_")
  
  plot.title = paste(tx, pheno.col, sep = " ")
  
  strains = sub("/", "_", hs.colors[,2])
  
  load(file = paste0(sanger.dir, "X.Rdata"))
  
  null.mod = glm(pheno[,pheno.col] ~ addcovar, family = gaussian)
  #null.mod = glm(trait ~ addcovar, family = poisson(link = "log"))
  null.ll = logLik(null.mod)
  pv = rep(0, nrow(sanger))
  
  glm.fxn = function(snp.rng, local.probs) {
    
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
    sdps.to.use = which(rowSums(cur.alleles) > 1.0)
    
    sex.col = which(colnames(addcovar) == "sex")
    if(length(sex.col) != 1) {
      stop("One of the columns of addcovar MUST be named 'sex'.")
    } # if(length(sex.col) != 1)
    
    # Run the model at each unique SDP.
    for(j in sdps.to.use) {
      
      
      full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = gaussian)
      #full.mod = glm(trait ~ addcovar + cur.alleles[j,], family = poisson(link = "log"))
      cur.ll[j] = logLik(full.mod)
      
    } # for(j)
    
    # This is the LRS.
    cur.ll = cur.ll - null.ll
    
    # Return the results.
    cur.ll[m]
    
  } # glm.fxn()
  
  # SNPs before the first marker.
  snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,1])
    
  } # if(length(snp.rng) > 0)
  
  # SNPs between Markers.
  for(i in 1:(nrow(obj$markers)-1)) {
    
    snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                      sanger.hdr$POS <= obj$markers[i+1,3])
    
    if(length(snp.rng) > 0) {
      
      # Take the mean of the haplotype probs at the surrounding markers.
      pv[snp.rng] = glm.fxn(snp.rng, (obj$probs[,,i] +
                                        obj$probs[,,i+1]) * 0.5)
      
    } # if(length(snp.rng) > 0)
    
  } # for(i)
  
  # SNPs after the last marker.
  snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])
    
  } # if(length(snp.rng) > 0)
  
  # Convert LRS to p-values using the chi-squared distribution.
  pv = pchisq(2 * pv, df = 1, lower.tail = FALSE)
  pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)
  
  maxPV = max(-log10(pv[,6]))
  print(paste0("Maximum LOD: ", maxPV), digits = 1)
  rm(maxPV)
  
  save(pv, file = paste0(file.prefix, "_chr", chr, ".Rdata"))
  
  png(paste0(file.prefix, "_chr", chr,".png"), width = 2600,
      height = 1200, res = 130)
  plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
  mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
  dev.off()
  
  # Return the positions and p-values.
  return(pv)
  rm(sanger, sanger.hdr)
  print(paste(round(difftime(Sys.time(), timechrx, units = 'hours'), digits = 2)))
  
} # GRSDgauss.xchr()
################################################################################
GRSDgauss.permsfast = function(obj, pheno, pheno.col, addcovar, tx, sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/") {
  chr = obj$markers[1,2]
  
  setwd(outdir)
  
  file.prefix = paste(tx, pheno.col, sep = "_")
  
  plot.title = paste(tx, pheno.col, sep = " ")
  
  strains = sub("/", "_", hs.colors[,2])
  
  load(file = paste0(sanger.dir, chr, ".Rdata"))
  
  null.mod = glm(pheno[,pheno.col] ~ addcovar, family = gaussian)
  #null.mod = glm(trait ~ addcovar, family = poisson(link = "log"))
  null.ll = logLik(null.mod)
  pv = rep(0, nrow(sanger))
  
  glm.fxn = function(snp.rng, local.probs) {
    
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
    sdps.to.use = which(rowSums(cur.alleles) > 1.0)
    # Run the model at each unique SDP.
    for(j in sdps.to.use) {
      
      
      full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = gaussian)
      #full.mod = glm(trait ~ addcovar + cur.alleles[j,], family = poisson(link = "log"))
      cur.ll[j] = logLik(full.mod)
      
    } # for(j)
    
    # This is the LRS.
    cur.ll = cur.ll - null.ll
    
    # Return the results.
    cur.ll[m]
    
  } # glm.fxn()
  
  # SNPs before the first marker.
  snp.rng = which(sanger.hdr$POS <= obj$markers[1,3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,1])
    
  } # if(length(snp.rng) > 0)
  
  # SNPs between Markers.
  for(i in 1:(nrow(obj$markers)-1)) {
    
    snp.rng = which(sanger.hdr$POS > obj$markers[i,3] &
                      sanger.hdr$POS <= obj$markers[i+1,3])
    
    if(length(snp.rng) > 0) {
      
      # Take the mean of the haplotype probs at the surrounding markers.
      pv[snp.rng] = glm.fxn(snp.rng, (obj$probs[,,i] +
                                        obj$probs[,,i+1]) * 0.5)
      
    } # if(length(snp.rng) > 0)
    
  } # for(i)
  
  # SNPs after the last marker.
  snp.rng = which(sanger.hdr$POS > obj$markers[nrow(obj$markers),3])
  if(length(snp.rng) > 0) {
    
    pv[snp.rng] = glm.fxn(snp.rng, obj$probs[,,nrow(obj$markers)])
    
  } # if(length(snp.rng) > 0)
  
  # Convert LRS to p-values using the chi-squared distribution.
  pv = pchisq(2 * pv, df = 1, lower.tail = FALSE)
  pv = data.frame(sanger.hdr, pv, stringsAsFactors = FALSE)
  
  
  
  # Return the positions and p-values.
  return(pv)
  
  
} # GRSDgauss.permsfast()
################################################################################
GRSDgauss.perms = function(perms, chr = 1:19, Xchr = FALSE,
                           pheno, pheno.col, probs, K, addcovar,
                           markers, snp.file, outdir = "~/Dropbox/Rstudio.cloud/WD",
                           tx = "", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/") {
  begin <- Sys.time()
  begin
  
  samples = intersect(rownames(pheno), rownames(probs))
  samples = intersect(samples, rownames(addcovar))
  samples = intersect(samples, rownames(K[[1]]))
  stopifnot(length(samples) > 0)
  
  pheno = pheno[samples,,drop = FALSE]
  addcovar = addcovar[samples,,drop = FALSE]
  probs = probs[samples,,,drop = FALSE]
  
  # DEFINE TRAIT #
  
  file.prefix = paste(tx, pheno.col, sep = "_")
  
  plot.title = paste(tx, pheno.col, sep = " ")
  print(paste(plot.title, "Permutation Analysis:", Sys.time()))
  
  trait = pheno[,pheno.col]
  
  # LOGISTIC REGRESSION MODEL #
  for(i in 1:length(K)) {
    K[[i]] = K[[i]][samples, samples]
  } # for(i)
  
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
  
  ##
  result = vector("list", length(data))
  names(result) = names(data)
  females = which(pheno$sex == "0")
  males = which(pheno$sex == "1")
  
  permutations = matrix(1, nrow = perms, ncol = 2, dimnames = list(1:perms, c("A", "X")))
  sanger.dir = sanger.dir
  for(p in 1:perms) {
    LODtime = Sys.time()
    print(p)
    
    
    new.order = rep(0, length(trait))
    new.order[females] = sample(females)
    new.order[males] = sample(males)
    
    log.perm = trait[new.order]
    pheno["new.col"] <- log.perm
    
    
    min.a.pv = 1
    
    for(i in 1:length(chr)) {
      result = GRSDgauss.permsfast(data[[i]], pheno = pheno, pheno.col = "new.col", addcovar, tx, sanger.dir)
      min.a.pv = min(min.a.pv, min(result$pv))
    } #for(i)
    
    min.a.pv = min(min.a.pv, min(result$pv))
    min.x.pv = 1
    
    # if(Xchr) {
    #   result = GRSDbinom.xchr.permsfast(data[["X"]], pheno = phenonew, pheno.col = "new.col", addcovar, tx, sanger.dir)
    #   min.x.pv = min(result$pv)
    # }
    
    # Save the minimum p-values.
    permutations[p,] = c(-log10(min.a.pv), -log10(min.x.pv))
    print(paste("Max random  LOD was", -log10(min.a.pv), -log10(min.x.pv)))
    print(paste(round(difftime(Sys.time(), LODtime, units = 'mins'), digits = 2),
                "minutes..."))
    
  }
  print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
              "hours elapsed during analysis"))
  
  save(permutations, file.prefix, file = paste0(file.prefix, "_perms.Rdata"))
  return(permutations)
  
  
} # GRSDgauss.perms

################################################################################
get.sig.thr = function(perms, alpha = 0.05, Xchr = F) {
  
  sig.thr = rep(0, length(alpha))
  
  if(Xchr) {
    
    if(!is.matrix(perms)) {
      stop(paste("'perms' is not a matrix. 'perms' must be a matrix",
                 "with 2 columns, named 'A' and 'X'."))
    } # if(!is.matrix(perms))
    
    if(!(all(colnames(perms) %in% c("A", "X")))) {
      stop(paste("The colnames of 'perms' are not equal to 'A' and",
                 "'X'. 'perms' must be a matrix, with 2 columns, named",
                 "'A' and 'X'."))
    } # if(!(all(colnames(perms) %in% c("A", "X"))))
    
    chrlen = get.chr.lengths()
    len.auto = sum(chrlen[1:19])
    len.X = chrlen["X"]
    len.all = len.auto + len.X
    alpha.auto = 1.0 - (1.0 - alpha)^(len.auto / len.all)
    alpha.X    = 1.0 - (1.0 - alpha)^(len.X / len.all)
    
    sig.thr = cbind("A" = quantile(perms[,"A"], probs = 1.0 - alpha.auto, na.rm = TRUE),
                    "X" = quantile(perms[,"X"], probs = 1.0 - alpha.X, na.rm = TRUE))
    rownames(sig.thr) = alpha
    
  } else {
    
    sig.thr = quantile(perms, probs = 1.0 - alpha, na.rm = TRUE)
    names(sig.thr) = alpha
    
  } # else
  
  return(sig.thr)
  
} # get.sig.thr()
################################################################################
