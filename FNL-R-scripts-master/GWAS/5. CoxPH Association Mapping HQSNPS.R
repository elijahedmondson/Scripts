# LOAD PACKAGES #

install.packages("C:/Users/edmondsonef/Desktop/DOQTL_1.0.5.tar.gz", repos = NULL, type ='source')

#'annotationTools', 'ggbio', 'hwriter', 'MUGAExampleData', 'QTLRel'

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("QTLRel")

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
outdir = "~/Desktop/files/"
options(stringsAsFactors = F)
setwd("~/Desktop/files/")
load(file = "~/Desktop/R/QTL/WD/GRSD.Rdata")
datadir <-"C:/Users/edmondsonef/Desktop/R-plots/"
setwd(datadir)

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
                   THB.merge = as.numeric(Total$Thyroid.HCC.Bmerge))
addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))

pheno["survival"] = rep(1, 1820)

#pheno = pheno[complete.cases(pheno$cat2),]

HZE = subset(pheno, group == "HZE")
Gamma = subset(pheno, group == "Gamma")
Unirradiated = subset(pheno, group == "Unirradiated")
All.irr = subset(pheno, Unirradiated == 0)


GRSD.coxph(pheno = HZE, pheno.col = "HCC.translocation", days.col = "days", probs, K, addcovar, markers, snp.file,
           outdir = "~/Desktop/files/", tx = "", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")


perms <- GRSDcoxph.perms(perms = 5, pheno = HZE, chr = 19, pheno.col = "Thyroid", days.col = "days",
                         probs = probs, K = K, addcovar = addcovar, markers = markers, snp.file = snp.file,
                         outdir = "~/Desktop/files/", tx = "HZE", sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/")




cox.bstrap = HS.cox.bootstrap(perms = 200, days.col = "days", tx = "", probs, K, addcovar, markers, snp.file, outdir = "~/Desktop/files/",
                              window = 8e+06, sanger.dir = "~/Desktop/R/QTL/WD/HS.sanger.files/",
                              chr = 4,
                              peakMB = 96173703,
                              pheno = Gamma,
                              pheno.col = "LSA.PreT")





# COX PH MODEL #
gamma = pheno[pheno$group == "Gamma",]
hze = pheno[pheno$group == "HZE",]
un = pheno[pheno$group == "Unirradiated",]

surv = Surv(un$days)
fit = survfit(surv ~ un$sex)
plot(fit, col = 1:2, las = 1, main = "plot.title")
legend("bottomleft", col = 1:2, lty = 1, legend = c("female", "male"))
mod = coxph(surv ~ un$sex)
text(x = 25, y = 0.15, labels = paste("p =", format(anova(mod)[2,4], digits = 2)), adj = 0)
fit


p1 <- ggplot(pheno, aes(x = family, y = days2))
p1 + geom_boxplot(notch = T, aes(fill = factor(family))) + geom_jitter() +
        theme_bw(base_size = 18) +
        theme(axis.text = element_text(size = 14),
              legend.position = "none",
              panel.grid.major = element_line(colour = "grey40"),
              panel.grid.minor = element_blank())


ggplot(cat4, aes(days4)) +
        geom_histogram(aes(y =..density.., fill = group), colour="black", binwidth = 14) +
        facet_grid(group ~ .) +
        theme_bw(base_size = 18) +
        geom_density(col=1) +
        ggtitle("Cataract Develop by Age\n(Cataract score 4.0)") +
        theme(legend.position = "none", axis.title.x = element_blank(), plot.margin=unit(c(1,1,1.5,1.2),"cm"))


multiplot(p1,p2, cols = 1)

p1 <- ggplot(pheno, aes(x=cat3, y=days3, colour=group)) +
        geom_line() +
        ggtitle("Growth curve for individual chicks")
p1

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
rm(Unirradiated, HZE, Gamma)

setwd(outdir)


