## qqman
## x = data.frame("BP," "CHR," "P," and optionally, "SNP")
library(qqman)
library(fastmatch)
library(Kmisc)
library(lattice)
library(HZE)

##### FOR QTL FROM WELLCOME TRUST ######

chr.function <- function(x){
        if(x == "chr1")
                return(1)
        if(x == "chr2")
                return(2)
        if(x == "chr3")
                return(3)
        if(x == "chr4")
                return(4)
        if(x == "chr5")
                return(5)
        if(x == "chr6")
                return(6)
        if(x == "chr7")
                return(7)
        if(x == "chr8")
                return(8)
        if(x == "chr9")
                return(9)
        if(x == "chr10")
                return(10)
        if(x == "chr11")
                return(11)
        if(x == "chr12")
                return(12)
        if(x == "chr13")
                return(13)
        if(x == "chr14")
                return(14)
        if(x == "chr15")
                return(15)
        if(x == "chr16")
                return(16)
        if(x == "chr17")
                return(17)
        if(x == "chr18")
                return(18)
        if(x == "chr19")
                return(19)
        else
                return(20)
}

qtl <- data.frame(SNP = gscan$scan$marker, CHR = gscan$scan$chromosome, BP = gscan$scan$bp, P = gscan$scan$logP)
qtl$CHR <- sapply(qtl$CHR, chr.function)
manhattan(qtl, main = 'HS: Context')

rm(gscan, qtl)


##### FOR QTL BINOM Function ######
plot.hs.qtl = function(qtl, bin.width = 1, ...) {

        new.qtl = NULL
        for(chr in 1:length(qtl)) {

                print(chr)

                # Create 100 SNP bins.
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
             col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
        mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)

        return(new.qtl)

} # plot.hs.qtl
load(file ="~/Desktop/R/QTL/WD/2.\ Binomial\ Mapping/Rdata/Allirr_Thyroid_QTL.Rdata")
Allirr = plot.hs.qtl(allirr)
load(file ="/Users/elijah/Desktop/R/QTL/WD/2.\ Binomial\ Mapping/Rdata/Gamma_Thyroid_GR_QTL.Rdata")
Gamma = plot.hs.qtl(gamma)
load(file ="/Users/elijah/Desktop/R/QTL/WD/2.\ Binomial\ Mapping/Rdata/HZE_Thyroid_GR_QTL.Rdata")
HZE = plot.hs.qtl(hze)


library(qqman)
HZE <- data.frame(CHR = HZE@seqnames,
                  BP = HZE@ranges@start,
                  P = HZE@elementMetadata@listData$p.value)

Gamma <- data.frame(CHR = Gamma@seqnames,
                      BP = Gamma@ranges@start,
                      P = Gamma@elementMetadata@listData$p.value)

All_Irradiated <- data.frame(CHR = Allirr@seqnames,
                      BP = Allirr@ranges@start,
                      P = Allirr@elementMetadata@listData$p.value)

HZE = HZE[which(HZE$CHR == 11), ]
Gamma = Gamma[which(Gamma$CHR == 11), ]
All_Irradiated = All_Irradiated[which(All_Irradiated$CHR == 11), ]

HZE$groups = c("HZE")
Gamma$groups = c("Gamma")
All_Irradiated$groups = c("All_Irradiated")

new = rbind(All_Irradiated, Gamma, HZE)
pv = new$P
bp = new$BP
chr = as.numeric(new$CHR)
groups = new$groups

manhattan_plot(pval = pv, bp, chr, groups, cex = 2.5, xlab = "", cutoff = c(5.73, 4.48))

manhattan_plot(pval = pv, bp, chr, groups, cex = 2, xlab = "")

chr.function <- function(x){
        if(x == "1")
                return(1)
        if(x == "2")
                return(2)
        if(x == "3")
                return(3)
        if(x == "4")
                return(4)
        if(x == "5")
                return(5)
        if(x == "6")
                return(6)
        if(x == "7")
                return(7)
        if(x == "8")
                return(8)
        if(x == "9")
                return(9)
        if(x == "10")
                return(10)
        if(x == "11")
                return(11)
        if(x == "12")
                return(12)
        if(x == "13")
                return(13)
        if(x == "14")
                return(14)
        if(x == "15")
                return(15)
        if(x == "16")
                return(16)
        if(x == "17")
                return(17)
        if(x == "18")
                return(18)
        if(x == "19")
                return(19)
        else
                return(20)
}
qtl$CHR <- sapply(qtl$CHR, chr.function)


qtl$CHR <- sapply(qtl$CHR, chr.function)





manhattan(qtl, main = 'Gamma Thyroid')
qq(qtl$P)

Allirr.cat2.13 <- Allirr.cat2[13]
Allirr.cat2.13 <- data.frame(BP = Allirr.cat2.13@unlistData@ranges@start, CHR = as.numeric(13), P = Allirr.cat2.13@unlistData@elementMetadata$p.value)
Gamma.cat2.13 <- Gamma.cat2[13]
Gamma.cat2.13 <- data.frame(BP = Gamma.cat2.13@unlistData@ranges@start, CHR = as.numeric(13), P = Gamma.cat2.13@unlistData@elementMetadata$p.value)
HZE.cat2.13 <- HZE.cat2[13]
HZE.cat2.13 <- data.frame(BP = HZE.cat2.13@unlistData@ranges@start, CHR = as.numeric(13), P = HZE.cat2.13@unlistData@elementMetadata$p.value)
UN.cat2.13 <- Unirradiated.cat2[13]
UN.cat2.13 <- data.frame(BP = UN.cat2.13@unlistData@ranges@start, CHR = as.numeric(13), P = UN.cat2.13@unlistData@elementMetadata$p.value)



manhattan(Allirr.cat2.13, ylim=c(0,10))
par(new=TRUE, bg = "transparent")
manhattan(Gamma.cat2.13, ylim=c(0,10))
par(new=TRUE, bg = "transparent")
manhattan(HZE.cat2.13, ylim=c(0,10))
manhattan(UN.cat2.13, ylim=c(0,10))

qq(Allirr.cat2.13$P)




##### FOR QTL FROM GRSDbinom ######

DOQTL:::plot.scanone.assoc(QTL.C1.train.frz, bin.size = 100, main = "HZE train_deltapctfrze_isi1_isi4: Cohort 1")

manhattan(new19, chr = "CHR", bp = "BP", p = "P",
           col = c("gray10", "gray60"), chrlabs = NULL,
           suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08),
           highlight = NULL, logp = TRUE)
qqman


# Look at a Quantile-Quantile plof of the p-values.
pv = unlist(sapply(qtl, function(z) { z$p.value }))
qqnorm(-log10(pv[sample(1:length(pv), 50000)]))
qqline(-log10(pv[sample(1:length(pv), 50000)]))
# This looks weird. Are there other covariates (like batch) that
# should be in your model? Or is the data full of 0 values with
# only a few non-zero values?


x <- Sys.time()
qq(new6$P)
Sys.time() - x

new1 <- data.frame(BP = qtl$`1`@ranges@start, CHR = as.numeric(1), P = qtl$`1`@elementMetadata$p.value)
new2 <- data.frame(BP = qtl$`2`@ranges@start, CHR = as.numeric(2), P = qtl$`2`@elementMetadata$p.value)
new3 <- data.frame(BP = qtl$`3`@ranges@start, CHR = as.numeric(3), P = qtl$`3`@elementMetadata$p.value)
new4 <- data.frame(BP = qtl$`4`@ranges@start, CHR = as.numeric(4), P = qtl$`4`@elementMetadata$p.value)
new5 <- data.frame(BP = qtl$`5`@ranges@start, CHR = as.numeric(5), P = qtl$`5`@elementMetadata$p.value)
new6 <- data.frame(BP = qtl$`6`@ranges@start, CHR = as.numeric(6), P = qtl$`6`@elementMetadata$p.value)
new7 <- data.frame(BP = qtl$`7`@ranges@start, CHR = as.numeric(7), P = qtl$`7`@elementMetadata$p.value)
new8 <- data.frame(BP = qtl$`8`@ranges@start, CHR = as.numeric(8), P = qtl$`8`@elementMetadata$p.value)
new9 <- data.frame(BP = qtl$`9`@ranges@start, CHR = as.numeric(9), P = qtl$`9`@elementMetadata$p.value)
new10 <- data.frame(BP = qtl$`10`@ranges@start, CHR = as.numeric(10), P = qtl$`10`@elementMetadata$p.value)
new11 <- data.frame(BP = qtl$`11`@ranges@start, CHR = as.numeric(11), P = qtl$`11`@elementMetadata$p.value)
new12 <- data.frame(BP = qtl$`12`@ranges@start, CHR = as.numeric(12), P = qtl$`12`@elementMetadata$p.value)
new13 <- data.frame(BP = qtl$`13`@ranges@start, CHR = as.numeric(13), P = qtl$`13`@elementMetadata$p.value)
new14 <- data.frame(BP = qtl$`14`@ranges@start, CHR = as.numeric(14), P = qtl$`14`@elementMetadata$p.value)
new15 <- data.frame(BP = qtl$`15`@ranges@start, CHR = as.numeric(15), P = qtl$`15`@elementMetadata$p.value)
new16 <- data.frame(BP = qtl$`16`@ranges@start, CHR = as.numeric(16), P = qtl$`16`@elementMetadata$p.value)
new17 <- data.frame(BP = qtl$`17`@ranges@start, CHR = as.numeric(17), P = qtl$`17`@elementMetadata$p.value)
new18 <- data.frame(BP = qtl$`18`@ranges@start, CHR = as.numeric(18), P = qtl$`18`@elementMetadata$p.value)
new19 <- data.frame(BP = qtl$`19`@ranges@start, CHR = as.numeric(19), P = qtl$`19`@elementMetadata$p.value)
new20 <- data.frame(BP = qtl$`X`@ranges@start, CHR = as.numeric(20), P = qtl$`X`@elementMetadata$p.value)

require(plyr)
x <- join_all(list(new1, new2, new3, new4, new5, new6, new7, new8, new9, new10, new11, new12, new13, new14,
           new15, new16, new17, new18, new19, new20), by = "BP", type = "full")
rm(new1, new2, new3, new4, new5, new6, new7, new8, new9, new10, new11, new12, new13, new14,
   new15, new16, new17, new18, new19, new20)
