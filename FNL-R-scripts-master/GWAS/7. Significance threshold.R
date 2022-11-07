
perms2 <- matrix(nrow = 1000, ncol = 2, -log10(perms[,1]), -log10(perms[,1]), dimnames = list(1:1000, c("A", "X")))

get.sig.thr(perms2)


perms = matrix(1, nrow = 1000, ncol = 2, dimnames = list(1:1000, c("A", "X")))

for(p in 1:1000) {
        print(p)
        pheno.perm = data.frame(row.names = sample(GRSD.pheno$row.names), sex = as.numeric(GRSD.pheno$sex == "M"),
                                cohort = as.numeric(GRSD.pheno$Cohort),
                                group = as.character(GRSD.pheno$groups),
                                family = as.character(GRSD.pheno$family),
                                weight = as.numeric(GRSD.pheno$Weight.corrected),
                                unirradiated = as.numeric(GRSD.pheno$Unirradiated),
                                Frz.total = as.numeric(GRSD.pheno$context_pctfrze_total),
                                pig.dis = as.numeric(GRSD.pheno$pigmentdispersion),
                                avgmo.total = as.numeric(GRSD.pheno$context_avgmo_total),
                                tone.frz = as.numeric(GRSD.pheno$cued_tone_pctfrze_total),
                                train.frz = as.numeric(GRSD.pheno$train_deltapctfrze_isi1_isi4),
                                train.shock = as.numeric(GRSD.pheno$train_deltaavgmot_shock1_shock5), 
                                Albino = as.numeric(GRSD.pheno$albino),
                                PulACA = as.numeric(GRSD.pheno$Pulmonary.Adenocarcinoma),
                                HCC = as.numeric(GRSD.pheno$Hepatocellular.Carcinoma),
                                LSA.PreT = as.numeric(GRSD.pheno$PreT))
        
        HZE.1perm <- subset(pheno.perm, group == "HZE" & cohort == 1)
        
        HZE.1addperm = matrix(HZE.1perm$sex, ncol = 1, dimnames = list(rownames(HZE.1perm), "sex"))
        
        min.a.pv = 1
        
        
        qtl = scanone.assoc(pheno = HZE.1perm, pheno.col = 11, probs = model.probs, K = K, addcovar = HZE.1addperm, markers = MM_snps, sdp.file = sdp.file, ncl = 4)
        
        min.a.pv = min(min.a.pv, min(qtl$`1`@elementMetadata$p.value), 
                       min(qtl$`2`@elementMetadata$p.value),
                       min(qtl$`3`@elementMetadata$p.value),
                       min(qtl$`4`@elementMetadata$p.value),
                       min(qtl$`5`@elementMetadata$p.value),
                       min(qtl$`6`@elementMetadata$p.value),
                       min(qtl$`7`@elementMetadata$p.value),
                       min(qtl$`8`@elementMetadata$p.value),
                       min(qtl$`9`@elementMetadata$p.value),
                       min(qtl$`10`@elementMetadata$p.value),
                       min(qtl$`11`@elementMetadata$p.value),
                       min(qtl$`13`@elementMetadata$p.value),
                       min(qtl$`14`@elementMetadata$p.value),
                       min(qtl$`15`@elementMetadata$p.value),
                       min(qtl$`16`@elementMetadata$p.value),
                       min(qtl$`17`@elementMetadata$p.value),
                       min(qtl$`18`@elementMetadata$p.value),
                       min(qtl$`19`@elementMetadata$p.value),
                       min(qtl$`X`@elementMetadata$p.value))
        
        print(-log10(min.a.pv))
        
        perms[p,] = c(-log10(min.a.pv), 1)
        
} # for(p)






#Significance Threshold Function
get.sig.thr = function(perms, alpha = 0.05, Xchr = TRUE) {
        
        sig.thr = rep(0, length(alpha))
        
        if(Xchr) {
                
                if(!is.matrix(perms)) {
                        stop(paste("'perms' is not a matrix. 'perms' must be a matrix",
                                   "with 2 columns, named 'A' and 'X'."))
                } # if(!is.matrix(perms))
                
                if(!(all(colnames(perms) %in% c("A", "X")))) {
                        stop(paste("The colnames of 'perms' are not equal to 'A' and",
                                   "'X'. 'perms' must be a matrix, with 2 columns, named",
                                   "'A' and 'X'."))
                } # if(!(all(colnames(perms) %in% c("A", "X"))))
                
                chrlen = get.chr.lengths()
                len.auto = sum(chrlen[1:19])
                len.X = chrlen["X"]
                len.all = len.auto + len.X
                alpha.auto = 1.0 - (1.0 - alpha)^(len.auto / len.all)
                alpha.X    = 1.0 - (1.0 - alpha)^(len.X / len.all)
                
                sig.thr = cbind("A" = quantile(perms[,"A"], probs = 1.0 - alpha.auto, na.rm = TRUE),
                                "X" = quantile(perms[,"X"], probs = 1.0 - alpha.X, na.rm = TRUE))
                rownames(sig.thr) = alpha
                
        } else {
                
                sig.thr = quantile(perms, probs = 1.0 - alpha, na.rm = TRUE)
                names(sig.thr) = alpha
                
        } # else
        
        return(sig.thr)
        
} # get.sig.thr()



