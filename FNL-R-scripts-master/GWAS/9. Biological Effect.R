
# LOAD PACKAGES #
library(BSgenome.Mmusculus.UCSC.mm10)
library(doParallel)
library(foreach)
library(Rsamtools)
library(VariantAnnotation)
library(DOQTL)
library(GenomicRanges)
library(survival)
library(regress)
library(HZE)
library(modEvA)
library(coxphw)
outdir = "~/Desktop/files/"
options(stringsAsFactors = F)
setwd("~/Desktop/files/")
load(file = "~/Desktop/R/QTL/WD/GRSD.Rdata")


# PHENOTYPE #
Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = Total$row.names, sex = as.numeric(Total$sex == "M"),
                   cohort = as.numeric(Total$Cohort),
                   Unirradiated = as.numeric(Total$Unirradiated),
                   family = as.numeric(Total$family),
                   group = as.character(Total$groups),
                   days = as.numeric(Total$days),
                   NN = as.numeric(Total$non.neoplastic),
                   Tumor = as.numeric(Total$Tumor),
                   albino = as.numeric(Total$albino),
                   black = as.numeric(Total$black),
                   grey = as.numeric(Total$grey),
                   brown = as.numeric(Total$creme.brown),
                   days2 = as.numeric(Total$Cataract.2.0.Score.Days),
                   cat2 = as.numeric(Total$Cataract.2.0.Score.Event),
                   days3 = as.numeric(Total$Cataract.3.0.Score.Days),
                   cat3 = as.numeric(Total$Cataract.3.0.Score.Event),
                   days4 = as.numeric(Total$Cataract.4.0.Score.Days),
                   cat4 = as.numeric(Total$Cataract.4.0.Score.Event),
                   pigdisp = as.numeric(Total$pigment.dispersion),
                   dilate = as.numeric(Total$Did.Not.Dilate),
                   PulACA = as.numeric(Total$Pulmonary.Adenocarcinoma),
                   HCC = as.numeric(Total$Hepatocellular.Carcinoma),
                   HSA = as.numeric(Total$Hemangiosarcoma),
                   HS = as.numeric(Total$Histiocytic.Sarcoma),
                   MammACA = as.numeric(Total$Mammary.Gland.Adenocarcinoma),
                   GCT = as.numeric(Total$Granulosa.Cell.Tumor),
                   Thyroid = as.numeric(Total$Thyroid.Tumor),
                   ThyroidAD = as.numeric(Total$Thyroid.Adenoma),
                   STS = as.numeric(Total$Soft.Tissue.Sarcomas),
                   AML = as.numeric(Total$Myeloid.Leukemia),
                   HardACA = as.numeric(Total$Harderian.Gland.Adenocarcinoma),
                   Harderian = as.numeric(Total$Harderian.Tumor),
                   HardAD = as.numeric(Total$Harderian.Gland.Adenoma),
                   LSA.BLL= as.numeric(Total$BLL),
                   LSA.Bmerge= as.numeric(Total$B.merge),
                   LSA.DLBCL= as.numeric(Total$DLBCL),
                   LSA.FBL= as.numeric(Total$FBL),
                   LSA.PreT = as.numeric(Total$PreT),
                   OSA = as.numeric(Total$Osteosarcoma),
                   PitAd = as.numeric(Total$Pituitary.Adenoma),
                   Amyloid = as.numeric(Total$Amyloidosis))
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))
pheno["survival"] = rep(1, 1820)
HZE = subset(pheno, group == "HZE")
Gamma = subset(pheno, group == "Gamma")
Unirradiated = subset(pheno, group == "Unirradiated")
All.irr = subset(pheno, Unirradiated == 0)




options(error=recover)
options(error=utils::recover) 
reach_full_in <- reachability(krack_full, 'in')



get.effect.size.coxPH(pheno = All.irr, pheno.col = "cat2", days.col = "days2", probs = probs, 
                      sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz", markers, 
                      threshold = 5.73, dir = "")






get.effect.size(pheno = Unirradiated, pheno.col = "cat2", probs = probs, sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz",
                markers, threshold = 5.73, dir = "~/Desktop/R/QTL/WD/7.\ Cataract/Logistic\ 2.0/BE/")




#### Alternative SNP method ####

get.effect.size2 = function(pheno = All.irr, pheno.col, probs = probs, sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz",
                           markers, threshold = 5.05, dir = "/Users/elijah/Desktop/R/QTL/WD/2.\ Binomial\ Mapping/")
{
        library(Rsamtools)

        #Enter the directory of QTL files
        files <- (Sys.glob(paste0(dir,"*.Rdata")))

        # Create a matrix of SDPs.
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

        for(j in 1:length(files)){
                #Create the matrix in which all data will be stored
                EFFECT = matrix(0, nrow = 20, ncol = 15,
                                dimnames = list(1:20, c("PHENOTYPE", "CHR", "SNP", "LOD", "ODDS",
                                                        "2.5% ODDS", "97.5% ODDS", "ANOVA Pr(>Chi)",
                                                        "AIC", "CoxSnell", "Nagelkerke", "McFadden",
                                                        "Tjur","sqPearson", "sqD")))
                load(file = files[j])
                print(files[j])
                for(i in 1:19) {

                        # Determine most significant SNP and LOD score on the chromosome of interest.
                        qtli = as.data.frame(qtl[[i]])
                        qtli = qtli[match(markers$Mb_NCBI38, qtli$start, nomatch=0),]
                        stopifnot(length(qtli) > 0)
                        SNP = qtli$start[which.min(qtli$p.value)]
                        LOD = -log10(qtli$p.value[which(qtli$start == SNP)])


                        # Run the loop for all significant loci
                        if(LOD > threshold) {
                                print(i)
                                # Read in the unique SDPs.
                                tf = TabixFile(file = sdp.file)
                                sdps = scanTabix(file = sdp.file, param = GRanges(seqnames = i, ranges = SNP))[[1]]
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

                                # Fit the model.
                                samples = intersect(rownames(pheno), rownames(probs))
                                samples = intersect(samples, rownames(addcovar))
                                samples = intersect(samples, rownames(geno))
                                stopifnot(length(samples) > 0)
                                pheno = pheno[samples,,drop = FALSE]
                                addcovar = addcovar[samples,,drop = FALSE]
                                geno = geno[samples,,drop = FALSE]
                                addcovar = addcovar[samples,,drop = FALSE]
                                #probs = probs[samples,,,drop = FALSE]

                                mod0 = glm(pheno[,pheno.col] ~ addcovar, family = binomial(logit))
                                #mod0
                                mod1 = glm(pheno[,pheno.col] ~ addcovar + geno[,1], family = binomial(logit))
                                #mod1
                                #summary(mod1)


                                ANOVA = anova(mod0,mod1,test = "Chisq")

                                #print("Odds of developing tumor as a function of genotype:")
                                odds = exp(coef(mod1))
                                #odds
                                oddsCI = exp(confint.default(mod1))
                                #oddsCI

                                #There are several ways of calculating (pseudo) R-squared values for logistic regression models,
                                #with no consensus about which is best. The RsqGLM function, now included in the modEvA package,
                                #calculates those of McFadden (1974), Cox & Snell (1989), Nagelkerke (1991), Tjur (2009), and
                                #the squared Pearson correlation between observed and predicted values.
                                R2 = RsqGLM(model = mod1)

                                #Linear models come with an R-squared value that measures the proportion of variation that the
                                #model accounts for. The R-squared is provided with summary(model) in R. For generalized linear
                                #models (GLMs), the equivalent is the amount of deviance accounted for D-squared (Guisan &
                                #Zimmermann 2000), but this value is not normally provided with the model summary. The Dsquared
                                #function, now included in the modEvA package (Barbosa et al. 2014), calculates it. There is also
                                #an option to calculate the adjusted D-squared, which takes into account the number of observations
                                #and the number of model parameters, thus allowing direct comparison among different models
                                #(Weisberg 1980, Guisan & Zimmermann 2000).
                                D2 = Dsquared(model = mod1)

                                EFFECT[i,] = c(pheno.col, chr, SNP, LOD, odds[3], oddsCI[3], oddsCI[6], ANOVA$`Pr(>Chi)`[2],
                                               mod1$aic, R2$CoxSnell, R2$Nagelkerke, R2$McFadden, R2$Tjur, R2$sqPearson, D2)
                                print(EFFECT[i,])
                                rm(LOD, oddsCI, ANOVA, mod1, mod0, R2, D2, SNP)
                        }

                }
                write.csv(EFFECT, file = paste0(files[j], "QTL", ".csv"))
                print(EFFECT)
                rm(qtl, EFFECT)
        }


}




get.effect.size.coxPH2 = function(pheno = pheno, pheno.col, days.col = days, probs = probs, sdp.file = "~/Desktop/R/QTL/WD/HS_Sanger_SDPs.txt.bgz",
                                 markers, threshold = 5.05, dir = "/Users/elijah/Desktop/R/QTL/WD/2.\ Binomial\ Mapping/")
{
        load("/Users/elijah/Desktop/R/QTL/WD/hs.colors.Rdata")
        #Enter the directory of QTL files
        files <- (Sys.glob(paste0(dir,"*.Rdata")))

        # Create a matrix of SDPs.
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

        for(j in 1:length(files)){
                #Create the matrix in which all data will be stored
                EFFECT = matrix(0, nrow = 20, ncol = 10,
                                dimnames = list(1:20, c("PHENOTYPE", "CHR", "SNP", "LOD", "Logrank P",
                                                        "BB Hazard", "BB Hazard Ratio", "SE BB HR", "Pr(>|z|)", "Rsquare")))
                load(file = files[j])
                print(files[j])
                for(i in 1:19) {
                        tryCatch({

                                # Determine most significant SNP and LOD score on the chromosome of interest.
                                qtli = as.data.frame(qtl[[i]])
                                qtli = qtli[match(markers$Mb_NCBI38, qtli$start, nomatch=0),]
                                stopifnot(length(qtli) > 0)
                                SNP = qtli$start[which.min(qtli$p.value)]
                                LOD = -log10(qtli$p.value[which(qtli$start == SNP)])


                                # Run the loop for all significant loci
                                if(LOD > threshold) {
                                        print(i)
                                        # Read in the unique SDPs.
                                        tf = TabixFile(file = sdp.file)
                                        sdps = scanTabix(file = sdp.file, param = GRanges(seqnames = i, ranges = SNP))[[1]]
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
                                        geno = round(geno, digits = 1)
                                        geno = ifelse(geno < 0.25, "AA",
                                                      ifelse(geno >=.25 & geno <= 0.75, "AB",
                                                             ifelse(geno > .75, "BB",
                                                                    NA)))


                                        # Fit the model.
                                        samples = intersect(rownames(pheno), rownames(probs))
                                        samples = intersect(samples, rownames(addcovar))
                                        samples = intersect(samples, rownames(geno))
                                        stopifnot(length(samples) > 0)
                                        pheno = pheno[samples,,drop = FALSE]
                                        addcovar = addcovar[samples,,drop = FALSE]
                                        geno = geno[samples,,drop = FALSE]
                                        addcovar = addcovar[samples,,drop = FALSE]
                                        #probs = probs[samples,,,drop = FALSE]

                                        ###### CoxPH model ######

                                        # Measures of explained variation, such as the coefficient of determination (R2) in linear models,
                                        # are helpful in assessing the explanatory power of a model. In survival analysis, these measures
                                        # help quantify the ability of prognostic factors to predict a patient's time until death. As in
                                        # linear models, covariates in Cox regression may be statistically significant but still have very
                                        # little predictive power. In the censored data setting, the definition of such a measure is not
                                        # straightforward; several measures of explained variation have been proposed. The most popular of
                                        # these is the generalized R-squared, calculated as 1-exp((χLR2)/n), where (χLR2) is the chi-square
                                        # statistic for the likelihood ratio test for the overall model, and n is the total number of
                                        # patients. Although the generalized R-squared is commonly recommended for the Cox model, its
                                        # sensitivity to the proportion of censored values is not often mentioned. In fact, the expected
                                        # value of R-squared decreases substantially as a function of the percent censored, with early
                                        # censoring having a greater impact than later censoring. Simulations show that complete data
                                        # R-squared values from the Cox model are very close to those from a similar linear model. However,
                                        # average R-squared values can decrease by 20% or more (e.g., R-squared from 0.5 to 0.4) with heavy
                                        # censoring (e.g., 50% censoring) compared to complete data. Simulation results will be presented,
                                        # and alternatives to the generalized R-squared will be discussed. SCHEMPER 1996

                                        surv = Surv(pheno[,days.col], pheno[,pheno.col])
                                        fit = survfit(surv ~ geno)
                                        mod = coxph(surv ~ geno)

                                        png(paste0(pheno.col, "_chr", chr,".png"), width = 2000,
                                            height = 1600, res = 250)
                                        plot(fit, col = 1:3, las = 1, main = paste0(pheno.col, ": chr ", chr, " bp ", SNP))
                                        legend("bottomleft", col = 1:3, lty = 1, legend = c("AA", "AB", "BB"))
                                        text(x = 5, y = 0.25, labels = paste("BB Hazard =", format(mod$coefficients[2], digits = 4),
                                                                             "(HR =", format(exp(mod$coefficients[2]), digits = 4), ")",
                                                                             "P.value =", format(anova(mod)[2,4], digits = 2)), adj = 0)
                                        dev.off()




                                        EFFECT[i,] = c(pheno.col, chr, SNP, LOD, format(anova(mod)[2,4]),
                                                       mod$coefficients[2], summary(mod)$coefficients["genoBB","exp(coef)"],
                                                       summary(mod)$coefficients["genoBB","se(coef)"],
                                                       summary(mod)$coefficients["genoBB","Pr(>|z|)"], summary(mod)$rsq[1])
                                        print(EFFECT[i,])
                                }
                        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})


                }
                write.csv(EFFECT, file = paste0(files[j], "QTL", ".csv"))
                print(EFFECT)
                rm(qtl, EFFECT)
        }


}#get.effect.size.coxPH()


