##### PLOTTING CoxPH and GLM GWAS #####
##### GAMMA == BLUE (#1F78B4(.122, .471, .706), #419fde(.255, .624, .871)) #####
##### HZE == RED (#DB2B3D (.859, .169, .239), #e66c79(.902, .424, .475)) #####
##### ALL.IRR ==  Purple (#663399 (.4, .2, .6), #8c53c6 (.549, .325, .776)) #####
##### UNIRRADIATED == GREEN (#33A02C (.2, .627, .173), #53ce4b (.325, .808 ,.294)) #####
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(DOQTL)
library(HZE)


par(mfrow = c(5,1), mar=c(1, 4, 1, 1) + 0.5)
#layout(matrix(4:1, 4, 1))
plot.hs.color.qtl(hze, bin.width = 100, color = "red", main = "HZE ion irradiated", ylim = c(0, 9.5))
abline(a = 5.73, b = 0, col = "grey")
#abline(a = 4.48, b = 0, col = "lightgrey")
plot.hs.color.qtl(gamma, bin.width = 100, color = "blue", main = "Gamma-ray irradiated", ylim = c(0, 9.5))
abline(a = 5.73, b = 0, col = "grey")
#abline(a = 4.48, b = 0, col = "lightgrey")
plot.hs.color.qtl(allirr, bin.width = 100, color = "purple", main = "All Irradiated", ylim = c(0, 9.5))
abline(a = 5.73, b = 0, col = "grey")
#abline(a = 4.48, b = 0, col = "lightgrey")
plot.hs.color.qtl(un, bin.width = 100, color = "green", main = "Unirradiated", ylim = c(0, 9.5))
abline(a = 5.73, b = 0, col = "grey")
#abline(a = 4.48, b = 0, col = "lightgrey")
plot.hs.color.qtl(all, bin.width = 100, color = "black", main = "All mice", ylim = c(0, 9.5))
abline(a = 5.73, b = 0, col = "grey")
#abline(a = 4.48, b = 0, col = "lightgrey")



par(mfrow = c(5,1), mar=c(1, 4, 1, 1) + 0.5)
#layout(matrix(4:1, 4, 1))
plot.hs.color.qtl(hze[15], bin.width = 10, color = "red", main = "HZE ion irradiated", ylim = c(0, 9.5))
abline(a = 5.73, b = 0, col = "grey")
#abline(a = 4.48, b = 0, col = "lightgrey")
plot.hs.color.qtl(gamma[15], bin.width = 10, color = "blue", main = "Gamma-ray irradiated", ylim = c(0, 9.5))
abline(a = 5.73, b = 0, col = "grey")
#abline(a = 4.48, b = 0, col = "lightgrey")
plot.hs.color.qtl(allirr[15], bin.width = 10, color = "purple", main = "All Irradiated", ylim = c(0, 9.5))
abline(a = 5.73, b = 0, col = "grey")
#abline(a = 4.48, b = 0, col = "lightgrey")
plot.hs.color.qtl(un[15], bin.width = 10, color = "green", main = "Unirradiated", ylim = c(0, 9.5))
abline(a = 5.73, b = 0, col = "grey")
#abline(a = 4.48, b = 0, col = "lightgrey")
plot.hs.color.qtl(all[15], bin.width = 10, color = "black", main = "All mice", ylim = c(0, 9.5))
abline(a = 5.73, b = 0, col = "grey")
#abline(a = 4.48, b = 0, col = "lightgrey")





plot.hs.color.qtl = function(qtl, bin.width = 1000, color = "black", ...) {

        new.qtl = NULL
        for(chr in 1:length(qtl)) {
                library(BSgenome.Mmusculus.UCSC.mm10)
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
        if(color == "purple") {
                col = rep(rgb(.4, .2, .6), length(new.qtl))
                even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
                col[even.chr] = rgb(.549, .325, .776)
                plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
                     col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
                mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)

        }
        if(color == "red") {
                col = rep(rgb(.859, .169, .239), length(new.qtl))
                even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
                col[even.chr] = rgb(.902, .424, .475)
                plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
                     col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
                mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)
        }
        if(color == "blue") {
                col = rep(rgb(.122, .471, .706), length(new.qtl))
                even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
                col[even.chr] = rgb(.255, .624, .871)
                plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
                     col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
                mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)
        }
        if(color == "green") {
                col = rep(rgb(.2, .627, .173), length(new.qtl))
                even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
                col[even.chr] = rgb(.325, .808 ,.294)
                plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
                     col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
                mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)

        }
        if(color == "black") {
                col = rep(rgb(0,0,0), length(new.qtl))
                even.chr = which(seqnames(new.qtl) %in% (1:10 * 2))
                col[even.chr] = rgb(0.7,0.7,0.7)
                plot(gmb, -log10(new.qtl$p.value), pch = 20, xaxt = "n",
                     col = col, las = 1, xlab = "", ylab = "-log10(p-value)", ...)
                mtext(side = 1, line = 0.5, at = chrmid, text = names(chrlen), cex = 1.2)

        }
        return(new.qtl)

} # plot.hs.color.qtl








library(qqman)
library(fastmatch)
library(Kmisc)
library(lattice)
library(HZE)



HZE <- data.frame(CHR = hze@unlistData@seqnames,
                  BP = hze@unlistData@ranges@start,
                  P = hze@unlistData@elementMetadata@listData$p.value)

Gamma <- data.frame(CHR = gamma@unlistData@seqnames,
                    BP = gamma@unlistData@ranges@start,
                    P = gamma@unlistData@elementMetadata@listData$p.value)

All_Irradiated <- data.frame(CHR = allirr@unlistData@seqnames,
                             BP = allirr@unlistData@ranges@start,
                             P = allirr@unlistData@elementMetadata@listData$p.value)

Unirradiated <- data.frame(CHR = un@unlistData@seqnames,
                           BP = un@unlistData@ranges@start,
                           P = un@unlistData@elementMetadata@listData$p.value)

HZE = HZE[which(HZE$CHR == 2), ]
Gamma = Gamma[which(Gamma$CHR == 2), ]
All_Irradiated = All_Irradiated[which(All_Irradiated$CHR == 2), ]
Unirradiated = Unirradiated[which(Unirradiated$CHR == 2), ]

HZE$groups = c("HZE")
Gamma$groups = c("Gamma")
All_Irradiated$groups = c("All_Irradiated")
Unirradiated$groups = c("Unirradiated")

new = rbind(All_Irradiated, HZE, Unirradiated, Gamma)
pv = new$P
bp = new$BP
chr = as.numeric(new$CHR)
groups = new$groups

manhattan_plot(pval = pv, bp, chr, groups, cex = 2.5, xlab = "", cutoff = c(5.73, 4.48))
