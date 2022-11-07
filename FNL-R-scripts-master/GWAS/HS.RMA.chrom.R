HS.RMA.chrom = function(perms, chr, pheno, pheno.col, probs, K, addcovar,
         markers, snp.file, outdir = "~/Desktop/",
         tx = "", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")
{
        begin <- Sys.time()
        file.prefix = paste0("Bootstrap.", tx, ".", pheno.col, "(Chr", chr, ")")
        plot.title = paste("Bootstrap ", tx, " ", pheno.col, " ", "(Chr ", chr, "):", sep = "")
        print(paste(plot.title, Sys.time()))
        
        load(file = paste0(sanger.dir, "probs", chr, ".Rdata"))
        
        samples = intersect(rownames(pheno), rownames(probs))
        samples = intersect(samples, rownames(addcovar))
        samples = intersect(samples, rownames(K[[1]]))
        stopifnot(length(samples) > 0)
        
        pheno = pheno[samples,,drop = FALSE]
        addcovar = addcovar[samples,,drop = FALSE]
        probs = probs[samples,,,drop = FALSE]
        
        #samples2 = markers$SNP_ID[which(markers$Chr == chr & markers$Mb_NCBI38 > (peakMB - 4000000) & markers$Mb_NCBI38 < (peakMB + 4000000))]
        samples2 = markers$SNP_ID[which(markers$Chr == chr)]
        probs = probs[,,samples2, drop = FALSE]
        
        K = K[[chr]][samples,samples,drop = FALSE]
        
        setwd(outdir)
        
        ##
        
        permutations = matrix(1, nrow = perms, ncol = 5, dimnames = list(1:perms, c("LOD", "min", "max", "average", "#Markers")))
        sanger.dir = sanger.dir
        rm(samples, samples2)
        for(p in 1:perms) {
                LODtime = Sys.time()
                print(p)
                
                #repeat {
                        phenoperm = data.frame(pheno[sample(nrow(pheno), replace = TRUE), ],
                                               check.names = FALSE, check.rows = FALSE)
                        
                      
                        
                        samples = sub("\\.[0-9]$", "", rownames(phenoperm))
                        probsperm = probs[samples,,]
                        rownames(probsperm) = make.unique(rownames(probsperm))
                        
                        
                        Kperm = K
                        Kperm = K[samples,samples,drop = FALSE]
                        rownames(Kperm) = make.unique(rownames(Kperm))
                        
                        
                        ### Move the model into the loop LOGISTIC REGRESSION MODEL ###
                        
                        Kperm = Kperm[samples, samples]
                        chrs = chr
                        data = vector("list", length(chrs))
                        names(data) = chrs
                        
                        rng = which(markers[,2] == chrs)
                        data = list(probsperm = probsperm, Kperm = Kperm,
                                    markers = markers[rng,])
                        
                        result = vector("list", length(data))
                        names(result) = names(data)
                        rm(probsperm, Kperm)
                        
                        
                        ### WORK HORSE ###
                        ### RETURN ANY POS FOR ENRIRE CHROMOSOME ###
                        
                        
                        result = GRSDbinom.permsfast(data, pheno = phenoperm, pheno.col, addcovar, tx, sanger.dir)
                        
                        top = max(-log10(result$pv))
                        MAX.LOD = result$POS[which(-log10(result$pv) == top)]
                        MegaBase = (min(MAX.LOD) + max(MAX.LOD))/2000000
                        
                        print(paste(MegaBase, "Mb: LOD", round(top, digits = 2)))
                        
                        #if (MAX.LOD > (peakMB - (window/2)) & MAX.LOD < (peakMB + (window/2))) break
                        
                #}
                
                
                
                print(paste0("Accepted locus: ", MegaBase, " Mb"))
                
                
                
                # Save the locus.
                permutations[p,] = c(top, min(MAX.LOD), max(MAX.LOD), ((min(MAX.LOD) + max(MAX.LOD))/2), length(MAX.LOD))
                print(paste(round(difftime(Sys.time(), LODtime, units = 'mins'), digits = 2),
                            "minutes..."))
                
        }
        
        print(paste(round(difftime(Sys.time(), begin, units = 'hours'), digits = 2),
                    "hours elapsed during analysis"))
        
        quant = quantile(permutations[,4], c(0.025,0.975))
        print(paste("95% Confidence Interval for QTL:", min(quant), "-", max(quant)))
        print(paste("Interval =", (round((((max(quant) - min(quant))/1000000)), digits = 2)), "Mb"))
        
        return(permutations)
        
        
}#HS.RMA.chrom()