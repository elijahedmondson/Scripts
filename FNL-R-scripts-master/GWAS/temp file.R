
plot.hs = function (x, chr, bin.size = 10000, ...)
{

        chrlen = get.chr.lengths()
        chrlen = chrlen[names(chrlen) %in% names(x)]
        chrlen = cumsum(chrlen)
        chrmid = c(0, chrlen[-length(chrlen)]) + diff(c(0, chrlen)) *
                0.5
        names(chrmid) = names(chrlen)
        pos = as.list(x@unlistData@ranges@start)
        #pos = lapply(x, FUN, start)
        pv = lapply(x, mcols)
        pv = lapply(pv, "[[", 1)
        for (c in 1:length(x)) {
                bins = seq(1, length(pv[[c]]), bin.size)
                bins[length(bins) + 1] = length(pv[[c]])
                pos2 = rep(0, length(bins))
                pv2 = rep(0, length(bins))
                for (i in 1:(length(bins) - 1)) {
                        wh = which.min(pv[[c]][bins[i]:bins[i + 1]])
                        wh = wh + bins[i] - 1
                        pos2[i] = pos[[c]][wh] * 1e-06 + max(1, chrlen[c -
                                                                               1])
                        pv2[i] = pv[[c]][wh]
                }
                pos[[c]] = pos2
                pv[[c]] = -log(pv2, 10)
        }
        chr = rep(names(x), sapply(pos, length))
        un.chr = unique(chr)
        m = match(chr, un.chr)
        col = m%%2 + 1
        plot(unlist(pos), unlist(pv), pch = 16, col = c("black",
                                                        "grey60")[col], las = 1, xaxt = "n", xlab = "", ylab = "-log10(p-value)",
             xaxs = "i")
        mtext(text = names(chrmid), side = 1, line = 0.1, at = chrmid)
}

