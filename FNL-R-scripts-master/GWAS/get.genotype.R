get.genotype = function(chr, pos, snp, markers, probs) {
        
        # See if the position is in bp or Mb.
        if(pos < 200) {
                pos = pos * 1e6
        } # if(pos < 200)
        
        # Convert the SNP to numbers.
        snp = unlist(snp)
        names(snp) = make.names(sub("_", ".", names(snp)))
        strains = make.names(hs.colors[,2])
        snp = snp[strains]
        snp = as.numeric(factor(snp)) - 1
        
        # Get the slices from the haplotype probs matrix.
        markers = markers[markers[,1] %in% dimnames(probs)[[3]],]
        probs = probs[,,dimnames(probs)[[3]] %in% markers[,1]]
        markers = markers[markers[,2] == chr,]
        probs = probs[,,markers[,1]]
        markers = markers[max(which(markers[,3] < pos * 1e-6)):min(which(markers[,3] > pos * 1e-6)),]
        
        # Get the probs for these markers.
        probs = probs[,,markers[,1], drop = FALSE]
        probs = apply(probs, 1:2, mean)
        
        # Multiply the two matrices and return the result.
        return(probs %*% snp)
        
} # get.genotype()

get.genotype(chr = 2, pos = 121483195, snp, markers, probs)

max.geno = get.max.geno(probs)


get.snp.details = function(gr, strains = character(), snp.file = 
                                   "ftp://ftp.jax.org/sanger/current_snps/mgp.v4.snps.dbSNP.vcf.gz") {
        
        hdr = scanVcfHeader(snp.file)
        
        if(length(strains) == 0) {
                strains = samples(hdr)
        } # if(length(strains) == 0)
        
        param = ScanVcfParam(info = "CSQ", geno = c("GT", "FI"), 
                             samples = strains, which = gr)
        
        snps = readVcf(file = snp.file, genome = "mm10", param = param)
        rowRanges(snps)$paramRangeID = names(rowRanges(snps))
        snps = genotypeCodesToNucleotides(snps)
        snps = expand(snps)
        
        alt = sapply(fixed(snps)$ALT, as.character)
        
        # Try to convert the consequences into one line.
        csq = info(snps)$CSQ
        keep = rep(FALSE, nrow(snps))
        for(i in 1:nrow(snps)) {
                
                spl = strsplit(csq[[i]], split = "\\|")
                spl = matrix(unlist(spl), ncol = length(spl))
                spl = spl[,spl[1,] == alt[i], drop = FALSE]
                
                tmp = paste0(spl[2,], spl[5,])
                spl = spl[,!duplicated(tmp), drop = FALSE]
                csq[[i]] = apply(spl, 2, paste, collapse = "|")
                
                keep[i] = alt[i] %in% unique(substring(sub("/", "", geno(snps)$GT[i,]), 1, 1))
                
        } # for(i)
        
        snps = snps[keep]
        
        snps
        
} # get.snp.details()


