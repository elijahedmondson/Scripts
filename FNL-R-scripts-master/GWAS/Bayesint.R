bayesint = function (qtl, chr, prob = 0.95, expandtomarkers = FALSE) 
{
        if (missing(qtl)) {
                stop("bayesint: The qtl cannot be null. Please supply a QTL object.")
        }
        if (missing(chr)) {
                stop(paste("bayesint: The chromosome cannot be null."))
        }
        else if (!chr %in% c(1:19, "X")) {
                stop(paste("bayesint: The chromosome must be 1 to 19 or X."))
        }
        if (prob < 0 | prob > 1) {
                stop(paste("bayesint: The probability must between 0 and 1."))
        }
        old.warn = options("warn")$warn
        options(warn = -1)
        if (is.numeric(as.numeric(chr))) {
                qtl = qtl$lod$A
        }
        else {
                qtl = qtl$lod$X
        }
        options(warn = old.warn)
        qtl[, 1] = as.character(qtl[, 1])
        qtl[, 2] = as.character(qtl[, 2])
        qtl[, 3] = as.numeric(qtl[, 3])
        qtl[, 7] = as.numeric(qtl[, 7])
        qtl = qtl[qtl[, 2] == chr, ]
        pos = qtl[, 3]
        if (any(is.na(pos))) {
                remove = which(is.na(pos))
                qtl = qtl[-remove, ]
                pos = pos[-remove]
        }
        breaks = approx(x = pos, y = 10^qtl[, 7], 
                        xout = seq(pos[1], pos[length(pos)], length.out = 1e+05))
        widths = diff(breaks$x)
        heights = breaks$y[-1] + breaks$y[-length(breaks$y)]
        trapezoids = 0.5 * heights * widths
        trapezoids = trapezoids/sum(trapezoids)
        ord = order(breaks$y[-length(breaks$y)], decreasing = TRUE)
        wh = min(which(cumsum(trapezoids[ord]) >= prob))
        int = range(ord[1:wh])
        left.snp = c(NA, qtl[1, 2], breaks$x[int][1], 
                     approx(qtl[, 3], qtl[, 4], breaks$x[int][1])$y, 
                     approx(qtl[, 3], qtl[, 5], breaks$x[int][1])$y, 
                     approx(qtl[, 3], qtl[, 6], breaks$x[int][1])$y, 
                     approx(qtl[, 3], qtl[, 7], breaks$x[int][1])$y)
        max.snp = qtl[which.max(qtl[, 7]), ]
        right.snp = c(NA, qtl[1, 2], breaks$x[int][2], 
                      approx(qtl[, 3], qtl[, 4], breaks$x[int][2])$y, 
                      approx(qtl[, 3], qtl[, 5], breaks$x[int][2])$y, 
                      approx(qtl[, 3], qtl[, 6], breaks$x[int][2])$y, 
                      approx(qtl[, 3], qtl[, 7], breaks$x[int][2])$y)
        if (expandtomarkers) {
                left.snp = qtl[max(which(breaks$x[int][1] >= qtl[, 3])), ]
                max.snp = qtl[which.max(qtl[, 7]), ]
                right.snp = qtl[min(which(breaks$x[int][2] <= qtl[, 3])), ]
        }
        retval = rbind(left.snp, max.snp, right.snp)
        retval[, 3] = round(as.numeric(retval[, 3]), digits = 6)
        retval[, 4] = round(as.numeric(retval[, 4]), digits = 6)
        retval[, 5] = round(as.numeric(retval[, 5]), digits = 6)
        retval[, 6] = round(as.numeric(retval[, 6]), digits = 6)
        retval$lod = as.numeric(retval[, 7])
        return(retval)
}

bayesint(qtl = qtl, chr = 2)
