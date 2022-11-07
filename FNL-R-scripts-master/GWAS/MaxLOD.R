
define.QTL(dir = "~/Desktop/files/", threshold = 5.05)

define.QTL = function(dir = "~/Desktop/files/completed/", threshold = 5.05) {
        
        files <- (Sys.glob(paste0(dir,"*.Rdata")))
        
        for(j in 1:length(files)){
                LODmat = matrix(0, nrow = 20, ncol = 5, dimnames = list(1:20, c("LOD", "peak", "1.5 min", "1.5 max", "interval range")))
                
                load(file = files[j])
                print(files[j])
                
                for(i in 1:19) {
                        tryCatch({
                                top = max(-log10(result[[i]]$pv))
                                if(top > threshold) {
                                        
                                        max.LOD.peak = result[[i]]$POS[which(-log10(result[[i]]$pv) == top)]
                                        
                                        LOD.drop.int = top - 1.5
                                        max.LOD.position = result[[i]]$POS[which(-log10(result[[i]]$pv) > max(LOD.drop.int))]
                                        
                                        print(paste(i, max.LOD.position[1], max(max.LOD.position), top))
                                        
                                        LODmat[i,] = c(top, ((min(max.LOD.peak) + max(max.LOD.peak))/2), 
                                                       max.LOD.position[1], max(max.LOD.position), 
                                                       (max(max.LOD.position) - max.LOD.position[1]))
                                }
                                }, error=function(e){})
                        
                        
                } 
                write.csv(LODmat, file = paste0(files[j], "QTL", ".csv"))
                rm(result)
        }
}




chr = 1

#Max LOD score
top <- max(-log10(qtl[[chr]]$p.value))
top
LOD.drop.int = top - 1.5
max.LOD.position <- qtl[[chr]]@ranges[which(-log10(qtl[[chr]]$p.value) == top)]
max.LOD.position

#Position of Max LOD
max.LOD.position <- qtl[[chr]]@ranges[which(-log10(qtl[[chr]]$p.value) == top)]
max.LOD.position

max.LOD.position <- qtl[[chr]]@ranges[which(-log10(qtl[[chr]]$p.value) > LOD.drop.int)]
max.LOD.position

start = max.LOD.position@start[1]
end = max(max.LOD.position@start[])

mgi = get.mgi.features(chr = chr, start = start, end = end, type = "gene", source = "MGI")
print(mgi$Name)



layout(matrix(3:1, 3, 1))
par(mfrow = c(2,2), mar=c(1, 4, 1, 1) + 0.1)
DOQTL:::plot.scanone.assoc(allirr, chr=4, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Gamma.days, chr=17, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Unirradiated.days, chr=17, bin.size = 100, main = "Unirradiated", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(qtl, chr=17, bin.size = 100, main = "Total Cases", ylim=c(0,15))

par(mfrow = c(3,1), mar=c(1, 4, 1, 1) + 0.5)
DOQTL:::plot.scanone.assoc(HZE.days, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Gamma.days, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Unirradiated.days, bin.size = 100, main = "Unirradiated", ylim=c(0,15))

par(mfrow = c(3,1), mar=c(1, 4, 1, 1) + 0.5)
DOQTL:::plot.scanone.assoc(HZE.days, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
abline(a = 13, b = 0, col = "red")
DOQTL:::plot.scanone.assoc(Gamma.days, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
abline(a = 13, b = 0, col = "red")
DOQTL:::plot.scanone.assoc(Unirradiated.days, bin.size = 100, main = "Unirradiated", ylim=c(0,15))
abline(a = 13, b = 0, col = "red")

qtl$coef$A[abs(qtl$coef$A[,7]) > 1,7]  = NA

coefplot = function (doqtl, chr = 1, stat.name = "LOD", conf.int = TRUE, 
          legend = TRUE, colors = "HS", sex, ...) 
{
        old.par = par(no.readonly = TRUE)
        cross = attr(doqtl, "cross")
        if (is.null(cross)) {
                if (colors[1] == "DO") {
                        colors = do.colors
                }
                else if (colors[1] == "HS") {
                        colors = hs.colors
                }
        }
        else {
                if (cross == "DO") {
                        colors = do.colors
                }
                else if (cross == "HS") {
                        colors = hs.colors
                }
        }
        num.founders = nrow(colors)
        call = match.call()
        lod = NULL
        coef = NULL
        if (chr == "X") {
                if (missing(sex)) {
                        stop("Sex (either M or F) must be specified on X chromosome.")
                }
                lod = doqtl$lod$X
                coef = doqtl$coef$X
                if (sex == "F") {
                        columns = match(paste("F", colors[, 1], sep = "."), 
                                        colnames(coef))
                }
                else {
                        columns = match(paste("M", colors[, 1], sep = "."), 
                                        colnames(coef))
                }
                columns = columns[!is.na(columns)]
                coef = coef[, c(1, columns)]
                colnames(coef)[1] = "A"
                colnames(coef) = sub("^[MF]\\.", "", colnames(coef))
                coef = coef[rownames(coef) %in% lod[, 1], ]
                coef[, 2:ncol(coef)] = coef[, 2:ncol(coef)] + coef[, 1]
                coef = coef - rowMeans(coef)
        }
        else {
                lod = doqtl$lod$A
                lod = lod[lod[, 2] == chr, ]
                intercept = doqtl$coef$A[, 1]
                coef = doqtl$coef$A[, (ncol(doqtl$coef$A) - num.founders + 
                                               1):ncol(doqtl$coef$A)]
                coef[, 1] = intercept
                colnames(coef)[1] = "A"
                coef = coef[rownames(coef) %in% lod[, 1], ]
                coef[, 2:ncol(coef)] = coef[, 2:ncol(coef)] + coef[, 1]
                coef = coef - rowMeans(coef)
        }
        if (!all(lod[, 1] == rownames(coef))) {
                stop(paste("The SNP IDs in column 1 of the qtl data frame must match", 
                           "the SNP IDs in the rownames of the coef matrix."))
        }
        if (!all(colnames(coef) %in% colors[, 1])) {
                stop(paste("The founder names in the colnames of the coefficient matrix", 
                           "must be in column 1 of the colors matrix."))
        }
        if (max(lod[, 3], na.rm = TRUE) > 200) {
                lod[, 3] = lod[, 3] * 1e-06
        }
        layout(mat = matrix(1:2, 2, 1), heights = c(0.66666667, 0.3333333))
        par(font = 2, font.lab = 2, font.axis = 2, las = 1, 
            plt = c(0.12, 0.99, 0, 0.85), xaxs = "i", lwd = 2)
        plot(lod[, 3], coef[, colors[1, 1]], type = "l", 
             col = colors[1, 3], lwd = 2, ylim = c(min(coef, na.rm = TRUE), max(coef * 
                                                                                                                                 2, na.rm = TRUE)), xlab = paste("Chr", chr), ylab = "Founder Effects", 
             axes = FALSE, ...)
        abline(v = 0:20 * 10, col = "grey80")
        for (i in 1:nrow(colors)) {
                points(lod[, 3], coef[, colors[i, 1]], type = "l", 
                       col = colors[i, 3], lwd = 2)
        }
        if (legend) {
                legend.side = "topleft"
                if (which.max(lod[, 7]) < nrow(lod)/2) {
                        legend.side = "topright"
                }
                legend(legend.side, colors[, 2], col = colors[, 3], lty = 1, 
                       lwd = 2, x.intersp = 0.75, y.intersp = 0.75, bg = "white", 
                       cex = 0.8)
        }
        axis(2)
        par(xpd = NA)
        usr = par("usr")
        rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
        par(xpd = FALSE)
        par(plt = c(0.12, 0.99, 0.3, 1))
        plot(lod[, 3], lod[, 7], type = "l", lwd = 2, xlab = "", 
             ylab = stat.name, ...)
        abline(v = 0:20 * 10, col = "grey80")
        points(lod[, 3], lod[, 7], type = "l", lwd = 2)
        if (conf.int) {
                interval = bayesint(doqtl, chr = chr)
                usr = par("usr")
                rect(interval[1, 3], usr[3], interval[3, 3], usr[4], 
                     col = rgb(0, 0, 1, 0.1), border = NA)
        }
        mtext(paste("Chr", chr), 1, 2)
        usr = par("usr")
        rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)
        par(old.par)
}


