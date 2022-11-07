#Regress a Peak



GRSD.assocGENO = function(pheno, pheno.col, chr, snp, probs, K, addcovar, markers, snp.file,
                      outdir = "~/Desktop/files/", tx = c("Gamma", "HZE", "Unirradiated", "All"),
                      sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/", sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz"){
        begin <- Sys.time()
        begin
        # COVARIATES #
        load("/Users/elijah/Desktop/R/QTL/WD/hs.colors.Rdata")
        samples = intersect(rownames(pheno), rownames(probs))
        samples = intersect(samples, rownames(addcovar))
        samples = intersect(samples, rownames(K[[1]]))
        stopifnot(length(samples) > 0)

        print(paste("A total of", length(samples), tx, "samples are complete."))

        pheno = pheno[samples,,drop = FALSE]
        addcovar = addcovar[samples,,drop = FALSE]
        probs = probs[samples,,,drop = FALSE]


        sdp.mat = matrix(as.numeric(intToBits(1:2^8)), nrow = 32)
        sdp.mat = sdp.mat[8:1,]
        dimnames(sdp.mat) = list(LETTERS[1:8], 1:2^8)

        #helper function from DG
        get.genotype = function(chr, pos, snp, markers, probs) {

                # Convert the SNP to numbers.
                snp = unlist(snp)
                names(snp) = make.names(sub("_", ".", names(snp)))
                strains = make.names(hs.colors[,2])

                # Get the slices from the haplotype probs matrix.
                markers = markers[markers[,1] %in% dimnames(probs)[[3]],]
                probs = probs[,,dimnames(probs)[[3]] %in% markers[,1]]
                markers = markers[markers[,2] == chr,]
                probs = probs[,,markers[,1]]
                markers = markers[max(which(markers[,3] < pos)):min(which(markers[,3] > pos)),]

                # Get the probs for these markers.
                probs = probs[,,markers[,1], drop = FALSE]
                probs = apply(probs, 1:2, mean)

                # Multiply the two matrices and return the result.
                return(probs %*% snp)

        } # get.genotype()


        # Read in the unique SDPs.
        tf = TabixFile(file = sdp.file)
        sdps = scanTabix(file = sdp.file, param = GRanges(seqnames = chr, ranges = snp))[[1]]
        sdps = strsplit(sdps, split = "\t")
        sdps = matrix(unlist(sdps), ncol = 3, byrow = T)
        chr  = sdps[1,1]
        pos  = as.numeric(sdps[,2])
        sdps = as.numeric(sdps[,3])

        geno = get.genotype(chr = chr,
                            pos = pos,
                            snp = sdp.mat[,sdps],
                            markers = markers,
                            probs = probs)
        #geno = round(geno, digits = 1)
        #geno = ifelse(geno < 0.25, "AA",
        #              ifelse(geno >=.25 & geno <= 0.75, "AB",
        #                     ifelse(geno > .75, "BB",
        #                            NA)))


        samples = intersect(rownames(pheno), rownames(probs))
        samples = intersect(samples, rownames(addcovar))
        samples = intersect(samples, rownames(geno))
        stopifnot(length(samples) > 0)
        pheno = pheno[samples,,drop = FALSE]
        addcovar = addcovar[samples,,drop = FALSE]
        geno = geno[samples,,drop = FALSE]
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

        result[[chr]] = GRSDbinom.regressGENO(data[[chr]], pheno, pheno.col, addcovar, tx, geno, snp)

}



GRSD.assocGENO(pheno = Gamma, pheno.col = "LSA.PreT", chr = 4, snp = 90248453, probs, K, addcovar, markers, snp.file,
               outdir = "~/Desktop/files/", tx = "Gamma",
               sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/", sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz")


