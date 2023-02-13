library(readxl)
library(BSgenome.Mmusculus.UCSC.mm10)
library(doParallel)
library(foreach)
library(Rsamtools)
#library(DOQTL)
#Bring in each function individually#
library(VariantAnnotation)
library(GenomicRanges)
library(regress)
library(MASS)
library(lmtest)
library(HZE)
library(ggplot2)
library(tidyverse)
library(gapminder)
library(dplyr)

#install.packages("C:/Users/edmondsonef/Desktop/DOQTL_0.99.1.tar.gz", repos = NULL, type="source")

# 1. GENOTYPE #

load(file = "C:/Users/edmondsonef/Desktop/QTL/HZEproject.Rdata")
load(file = "C:/Users/edmondsonef/Desktop/QTL/HS.colors.Rdata")
# load(file = "~/Desktop/R/QTL/WD/HS\ HMM\ Rdata/model.probs.Rdata")
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))
outdir = "C:/Users/edmondsonef/Desktop/R-plots/"


sdp.file = "C:/Users/edmondsonef/Desktop/QTL/HS_Sanger_SDPs.txt.bgz"


# 2. PHENOTYPE #
Total <- read_excel("C:/Users/edmondsonef/Desktop/Cataract/CATARACT_final.xlsx")
#Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = Total$row.names, sex = as.numeric(Total$sex == "M"),
                   albino = as.numeric(Total$`coat color` == 'albino'),
                   family = as.numeric(Total$family),
                   group = as.character(Total$groups),
                   days = as.numeric(Total$days),
                   cat = as.numeric(Total$cat_score))
                   # unirradiated = as.numeric(Total$Unirradiated),
                   # AML.t = as.numeric(Total$AML.transform),
                   # HCC.t = as.numeric(Total$HCC.transform),
                   # PreT.t = as.numeric(Total$PreT.transform),
                   # LSA.t = as.numeric(Total$LSA.transform),
                   # Amyloid.t = as.numeric(Total$Amyloid.transform),
                   # PitAd.t = as.numeric(Total$PitAd.transform),
                   # PulACA.t = as.numeric(Total$PulACA.transform),
                   # Hard.number.t = as.numeric(Total$Hard.number.transform),
                   # Hard.t = as.numeric(Total$Hard.transform),
                   # FBL.t = as.numeric(Total$FBL.transform),
                   # Bmerge.t = as.numeric(Total$Merge.Transform),
                   # DLBCL.t = as.numeric(Total$DLBCL.transform),
                   # NN.t = as.numeric(Total$NN.transform),
                   # PulMet.t = as.numeric(Total$PulMet.transform),
                   # Thyroid.t = as.numeric(Total$Thyroid.transform),
                   # OSA.t = as.numeric(Total$OSA.transform),
                   # PSC = as.numeric(Total$Pulmonary.Sarcomatoid.Carcinoma),
                   # unirradiated = as.numeric(Total$Unirradiated),
                   # Frz.total = as.numeric(Total$context_pctfrze_total),
                   # pig.dis = as.numeric(Total$pigment.dispersion),
                   # avgmo.total = as.numeric(Total$context_avgmo_total),
                   # tone.frz = as.numeric(Total$cued_tone_pctfrze_total),
                   # train.frz = as.numeric(Total$train_deltapctfrze_isi1_isi4),
                   # train.shock = as.numeric(Total$train_deltaavgmot_shock1_shock5),
                   # pu.1 = as.numeric(Total$Pu.1.transform),
                   # no.pu.1 = as.numeric(Total$No.deletion.AML.transform))

HZE <- subset(pheno, group == "HZE")
Gamma <- subset(pheno, group == "Gamma")
Un <- subset(pheno, group == "Unirradiated")
Allirr <- subset(pheno, group != "Unirradiated")


# 3. COVARIATES #

addcovar <- matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))



# 4. ASSOCIATION MAPPING #
qtl <- GRSD.assoc(pheno = Gamma, pheno.col = 2, probs, K, addcovar, markers, snp.file,
                      outdir = "C:/Users/edmondsonef/Desktop/", tx = "Gamma",
                      sanger.dir = "C:/Users/edmondsonef/Desktop/R/QTL/HS.sanger.files/")


qtl <- scanone(pheno = Gamma, pheno.col = 2, probs = probs, K = K, addcovar = addcovar, snps = MM_snps)


qtl <- scanone.assoc(pheno = HZE, pheno.col = 2, probs = probs, K = K, cross = HS,
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
