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
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))

HZE <- subset(pheno, group == "HZE")
Gamma <- subset(pheno, group == "Gamma")
Un <- subset(pheno, group == "Unirradiated")
Allirr <- subset(pheno, group != "Unirradiated")


# 3. COVARIATES #

addcovar <- matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))



# 4. ASSOCIATION MAPPING #
qtl <- GRSD.assoc(pheno = Gamma, pheno.col = 3, probs, K, addcovar, markers, snp.file,
                  outdir = "C:/Users/edmondsonef/Desktop/R-plots/", tx = "Gamma",
                  sanger.dir = "C:/Users/edmondsonef/Desktop/QTL/HS.sanger.files/")


qtl <- scanone(pheno = Gamma, pheno.col = 3, probs = probs, K = K, addcovar = addcovar, snps = markers)


qtl <- scanone.assoc(pheno = HZE, pheno.col = 3, probs = probs, K = K, cross = HS,
                     addcovar = addcovar, markers = MM_snps, sdp.file = sdp.file, ncl = 4)

DOQTL:::plot.scanone.assoc(qtl, bin.size = 100, main = "")



perms <- Scanone.assoc.perms(perms = 200, pheno = Gamma, pheno.col = "AML.t", probs = model.probs, 
                             K = K, tx = "Gamma", addcovar = addcovar, markers = MM_snps, 
                             sdp.file = sdp.file, ncl = 4)







png("CorrectNeuro.albino.png", width = 2400, height = 1080, res = 200)
DOQTL:::plot.scanone.assoc(qtl, bin.size = 100)
dev.off()

png("_CHR_4.png", width = 2000, height = 1600, res = 128)
DOQTL:::plot.scanone.assoc(qtl, chr = 7, bin.size = 100)
dev.off()

layout(matrix(3:1, 3, 1))
par(mfrow = c(2,2), mar=c(1, 4, 1, 1) + 0.1)
loop.hs.qtl(HZE.cat2, chr=17, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Gamma.cat2, chr=17, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Unirradiated.cat2, chr=17, bin.size = 100, main = "Unirradiated", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Allirr.cat2, chr=17, bin.size = 100, main = "All irradiated", ylim=c(0,15))

par(mfrow = c(3,1), mar=c(1, 4, 1, 1) + 0.5)
DOQTL:::plot.scanone.assoc(HZE.days, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
abline(a = 13, b = 0, col = "red")
DOQTL:::plot.scanone.assoc(Gamma.days, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
abline(a = 13, b = 0, col = "red")
DOQTL:::plot.scanone.assoc(Unirradiated.days, bin.size = 100, main = "Unirradiated", ylim=c(0,15))
abline(a = 13, b = 0, col = "red")



# 5. LINKAGE MAPPING #

qtl = scanone(pheno = HZE.1, pheno.col = "albino", probs = probs, K = K,
              addcovar = addcovar, snps = MM_snps)
plot(qtl, main = "")

perms = scanone.perm(pheno = pheno, pheno.col = "AML", probs = model.probs, addcovar = addcovar,
                     snps = MM_snps, path = "~/Desktop/",
                     nperm = 100)

thr1 = quantile(perms, probs = 0.90)
thr2 = quantile(perms, probs = 0.95)
thr3 = quantile(perms, probs = 0.99)

plot(qtl, sig.thr = c(thr1, thr2, thr3), main = "PSC")

interval = bayesint(qtl, chr = 7)
interval
mgi = get.mgi.features(chr = interval[1,2], start = interval[1,3],
                       end = interval[3,3], type = "gene", source = "MGI")
nrow(mgi)
head(mgi)

ma = assoc.map(pheno = pheno, pheno.col = "albino", probs = probs, K = K, addcovar = addcovar,
               snps = MM_snps, chr = interval[1,2], start = interval[1,3], end = interval[3,3])
coefplot(qtl, chr = 7, cross = "HS", colors = "HS")
tmp = assoc.plot(ma, thr = 1)
unique(tmp$sdps)





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

layout(matrix(3:1, 3, 1))
par(mfrow = c(2,2), mar=c(1, 4, 1, 1) + 0.1)
loop.hs.qtl(HZE.cat2, chr=17, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
loop.hs.qtl(Gamma.cat2, chr=17, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
loop.hs.qtl(Unirradiated.cat2, chr=17, bin.size = 100, main = "Unirradiated", ylim=c(0,15))
loop.hs.qtl(Allirr.cat2, chr=17, bin.size = 100, main = "All irradiated", ylim=c(0,15))


# Plotting for scanone.assoc.
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
  print(paste0("Maximum LOD: ", maxPV))
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
