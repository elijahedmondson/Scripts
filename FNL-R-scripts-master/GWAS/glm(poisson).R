

GRSD.assoc1 = function(pheno, pheno.col, probs, K, addcovar, markers, snp.file,
         outdir = "~/Desktop/", tx = c("Gamma", "HZE", "Unirradiated")){
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
        
        for(i in 1:19) {
                print(paste("CHROMOSOME", i))
                result[[i]] = GRSDbinom1(data[[i]], pheno, pheno.col, addcovar, tx)
        } #for(i)
        
        print("X CHROMOSOME")
        result[["X"]] = GRSDbinom.xchr1(data[["X"]], pheno, pheno.col, addcovar, tx)
        
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
        
        save(qtl, file.prefix, file = paste(tx, pheno.col, sep = ".", "_QTL.Rdata"))
        
        
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
GRSDbinom1 = function(obj, pheno, pheno.col, addcovar, tx) {
        
        chr = obj$markers[1,2]
        
        setwd(outdir)
        
        file.prefix = paste(tx, pheno.col, sep = "_")
        
        plot.title = paste(tx, pheno.col, sep = " ")
        
        strains = sub("/", "_", hs.colors[,2])
        
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
        
        
        sanger.hdr = data.frame(ID = names(rowRanges(sanger)), CHR = as.character(seqnames(sanger)),
                                POS = start(sanger), REF = as.character(fixed(sanger)$REF),
                                ALT = alt, stringsAsFactors = FALSE)
        rm(alt)
        
        
        sanger = cbind(geno(sanger)$GT[,1:4,drop = FALSE],
                       "C57BL_6J" = "0/0",
                       geno(sanger)$GT[,5:7,drop = FALSE])
        
        sanger = (sanger != "0/0") * 1
        
        # Make the MAF between 1/8 and 4/8.
        flip = which(rowSums(sanger) > 4)
        sanger[flip,] = 1 - sanger[flip,,drop = FALSE]
        rm(flip)
        
        #null.mod = glm(pheno[,pheno.col] ~ addcovar, family = binomial(logit))
        null.mod = glm(pheno[,pheno.col] ~ addcovar, family = poisson())
        #null.mod = glm(pheno[,pheno.col] ~ addcovar, family = gaussian())
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
                        full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = poisson())
                        #full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = gaussian())
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
        
        png(paste0(file.prefix, "_chr", chr,".png"), width = 2000,
            height = 1600, res = 200)
        plot(as.numeric(pv[,3]) * 1e-6, -log10(pv[,6]), pch = 20)
        mtext(side = 3, line = 0.5, text = paste(plot.title, ": Chr", chr))
        dev.off()
        
        # Return the positions and p-values.
        return(pv)
        
}
GRSDbinom.xchr1 = function(obj, pheno, pheno.col, addcovar, tx) {
        
        chr = obj$markers[1,2]
        
        setwd(outdir)
        
        file.prefix = paste(tx, pheno.col, sep = "_")
        
        plot.title = paste(tx, pheno.col, sep = " ")
        
        strains = sub("/", "_", hs.colors[,2])
        
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
        
        
        sanger.hdr = data.frame(ID = names(rowRanges(sanger)), CHR = as.character(seqnames(sanger)),
                                POS = start(sanger), REF = as.character(fixed(sanger)$REF),
                                ALT = alt, stringsAsFactors = FALSE)
        rm(alt)
        
        
        sanger = cbind(geno(sanger)$GT[,1:4,drop = FALSE],
                       "C57BL_6J" = "0/0",
                       geno(sanger)$GT[,5:7,drop = FALSE])
        
        sanger = (sanger != "0/0") * 1
        
        # Make the MAF between 1/8 and 4/8.
        flip = which(rowSums(sanger) > 4)
        sanger[flip,] = 1 - sanger[flip,,drop = FALSE]
        rm(flip)
        
        #null.mod = glm(pheno[,pheno.col] ~ addcovar, family = binomial(logit))
        null.mod = glm(pheno[,pheno.col] ~ addcovar, family = poisson())
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
                        full.mod = glm(pheno[,pheno.col] ~ addcovar + cur.alleles[j,], family = poisson())
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
        
}

