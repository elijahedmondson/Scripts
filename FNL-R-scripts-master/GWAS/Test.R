#Loading Packages and Data
library(DOQTL)
library(doParallel)
library(foreach)
library(Rsamtools)
library(VariantAnnotation)
library(GenomicRanges)
library(survival)
library(regress)
sample <- read.csv(file = "/Users/elijahedmondson/Desktop/sample.csv")
load(file = "/Users/elijahedmondson/Desktop/GRSD_master.Rdata")
cross = "HS"
options(stringsAsFactors = F)
snp.file = "/Users/elijahedmondson/Desktop/R/Build/mgp.v4.snps.dbSNP.vcf.gz"
sdp.file = "/Users/elijahedmondson/Desktop/R/Build/HS_Sanger_SDPs.txt.bgz"
setwd("/Users/elijahedmondson/Desktop/R/QTL/WD")
pheno = data.frame(row.names = sample$row.names, sex = as.numeric(sample$sex == "M"),
                   albino = as.numeric(sample$albino),
                   days = as.numeric(sample$days))
covar = data.frame(sex = as.numeric(sample$sex == "M"))
addcovar = covar
rownames(covar) = rownames(pheno)
rownames(addcovar) = rownames(pheno)

#Linkage Mapping
LM.qtl = scanone(pheno = pheno, pheno.col = "days", probs = model.probs, K = K,
              addcovar = covar, snps = MM_snps)
plot(LM.qtl, main = "test")
save(LM.qtl, file = "")
perms = scanone.perm(pheno = pheno, pheno.col = "days", probs = model.probs, addcovar = covar,
                     snps = MM_snps, path = "/Users/elijahedmondson/Desktop/R/QTL/perms",
                     nperm = 10)
thr1 = quantile(perms, probs = 0.90)
thr2 = quantile(perms, probs = 0.95)
thr3 = quantile(perms, probs = 0.99)
plot(qtl, sig.thr = c(thr1, thr2, thr3), main = "")
interval = bayesint(qtl, chr = 7)
interval
mgi = get.mgi.features(chr = interval[1,2], start = interval[1,3], end = interval[3,3],
                       type = "gene", source = "MGI")
nrow(mgi)
head(mgi)

#Association Mapping
pheno[,2] = as.numeric(pheno[,2]) - 1

AM.qtl = scanone.assoc(pheno = pheno, pheno.col = 2, probs = probs, K = K,
                    addcovar = covar, markers = markers, cross = "HS", sdp.file = sdp.file, ncl = 2)
plot(AM.qtl, main = "test")
save(AM.qtl, file = "")

png("albino_QTL.png", width = 2000, height = 1600, res = 128)
DOQTL:::plot.scanone.assoc(qtl, bin.size = 100)
dev.off()

png("albino_QTL_chr7.png", width = 2000, height = 1600, res = 128)
DOQTL:::plot.scanone.assoc(qtl, chr = 7, bin.size = 10)
dev.off()

















stack = data.frame(Gamma = c("Non-neoplastic" = sum(Gamma$non.neoplastic),
                               "Lymphoma PreT" = sum(Gamma$PreT),
                               "Lymphoma BLL" = sum(Gamma$BLL),
                               "Lymphoma FBL" = sum(Gamma$FBL),
                               "Lymphoma DLBCL" = sum(Gamma$DLBCL),
                               "Myeloid Leukemia" = sum(Gamma$Myeloid.Leukemia),
                               "Pulmonary Adenocarcinoma" = sum(Gamma$Pulmonary.Adenocarcinoma),
                               "Hepatocellular Carcinoma" = sum(Gamma$Hepatocellular.Carcinoma),
                               "Hemangiosarcoma" = sum(Gamma$Hemangiosarcoma),
                               "Histiocytic Sarcoma" = sum(Gamma$Histiocytic.Sarcoma),
                               "Mammary Adenocarcinoma" = sum(Gamma$Mammary.Gland.Adenocarcinoma),
                               "Ovarian Granulosa Cell Tumor" = sum(Gamma$Granulosa.Cell.Tumor),
                               "Thyroid Adenoma" = sum(Gamma$Thyroid.Tumor),
                               "Soft Tissue Sarcoma" = sum(Gamma$Soft.Tissue.Sarcomas),
                               "Harderian Gland Adenocarcinoma" = sum(Gamma$Harderian.Gland.Adenocarcinoma),
                               "Harderian Gland Adenoma" = sum(Gamma$Harderian.Gland.Adenoma),
                               "Osteosarcoma" = sum(Gamma$Osteosarcoma),
                               "Pituitary Adenoma" = sum(Gamma$Pituitary.Adenoma)),
                     HZE = c("Non-neoplastic" = sum(HZE$non.neoplastic),
                           "Lymphoma PreT" = sum(HZE$PreT),
                          "Lymphoma BLL" = sum(HZE$BLL),
                          "Lymphoma FBL" = sum(HZE$FBL),
                          "Lymphoma DLBCL" = sum(HZE$DLBCL),
                          "Myeloid Leukemia" = sum(HZE$Myeloid.Leukemia),
                          "Pulmonary Adenocarcinoma" = sum(HZE$Pulmonary.Adenocarcinoma),
                          "Hepatocellular Carcinoma" = sum(HZE$Hepatocellular.Carcinoma),
                          "Hemangiosarcoma" = sum(HZE$Hemangiosarcoma),
                          "Histiocytic Sarcoma" = sum(HZE$Histiocytic.Sarcoma),
                          "Mammary Adenocarcinoma" = sum(HZE$Mammary.Gland.Adenocarcinoma),
                          "Ovarian Granulosa Cell Tumor" = sum(HZE$Granulosa.Cell.Tumor),
                          "Thyroid Adenoma" = sum(HZE$Thyroid.Tumor),
                          "Soft Tissue Sarcoma" = sum(HZE$Soft.Tissue.Sarcomas),
                          "Harderian Gland Adenocarcinoma" = sum(HZE$Harderian.Gland.Adenocarcinoma),
                          "Harderian Gland Adenoma" = sum(HZE$Harderian.Gland.Adenoma),
                          "Osteosarcoma" = sum(HZE$Osteosarcoma),
                          "Pituitary Adenoma" = sum(HZE$Pituitary.Adenoma)),
                   Unirradiated = c("Non-neoplastic" = sum(Unirradiated$non.neoplastic),
                                    "Lymphoma PreT" = sum(Unirradiated$PreT),
                                    "Lymphoma BLL" = sum(Unirradiated$BLL),
                                    "Lymphoma FBL" = sum(Unirradiated$FBL),
                                    "Lymphoma DLBCL" = sum(Unirradiated$DLBCL),
                                    "Myeloid Leukemia" = sum(Unirradiated$Myeloid.Leukemia),
                                    "Pulmonary Adenocarcinoma" = sum(Unirradiated$Pulmonary.Adenocarcinoma),
                                    "Hepatocellular Carcinoma" = sum(Unirradiated$Hepatocellular.Carcinoma),
                                    "Hemangiosarcoma" = sum(Unirradiated$Hemangiosarcoma),
                                    "Histiocytic Sarcoma" = sum(Unirradiated$Histiocytic.Sarcoma),
                                    "Mammary Adenocarcinoma" = sum(Unirradiated$Mammary.Gland.Adenocarcinoma),
                                    "Ovarian Granulosa Cell Tumor" = sum(Unirradiated$Granulosa.Cell.Tumor),
                                    "Thyroid Adenoma" = sum(Unirradiated$Thyroid.Tumor),
                                    "Soft Tissue Sarcoma" = sum(Unirradiated$Soft.Tissue.Sarcomas),
                                    "Harderian Gland Adenocarcinoma" = sum(Unirradiated$Harderian.Gland.Adenocarcinoma),
                                    "Harderian Gland Adenoma" = sum(Unirradiated$Harderian.Gland.Adenoma),
                                    "Osteosarcoma" = sum(Unirradiated$Osteosarcoma),
                                    "Pituitary Adenoma" = sum(Unirradiated$Pituitary.Adenoma)))

library(RColorBrewer)


my.col <- c("#ccf2ff", "#9CB071", "#9CB071",
         "#9CB071", "#9CB071", "#87AFC7",
         "#FFCBA4", "#ff9966", "#E6E600FF",
         "#2B547E", "#E8C034FF", "#404040",
         "#e6e6e6", "#ff6666", "#617C58",
         "#87CEEB", "#C48793", "#EDC9AF")
barplot(t(stack), col = my.col, ylab = "Number of Mice", ylim = c(0,800), xlim = c(0,12), width = 2)


legend("bottomright",
       legend = c(colnames(stack)[18:1]), #in order from top to bottom
       fill = my.col[18:1], # 6:1 reorders so legend order matches graph
       title = "Tumor Histotype")







Group <- c(rep(c("Gamma", "HZE", "Unirradiated"), each = 18))
Histotype <- c(rep(c("Non-neoplastic", "Lymphoma PreT","Lymphoma BLL",
               "Lymphoma FBL", "Lymphoma DLBCL", "Myeloid Leukemia",
               "Pulmonary Adenocarcinoma", "Hepatocellular Carcinoma",
               "Hemangiosarcoma", "Histiocytic Sarcoma",
               "Mammary Adenocarcinoma",
               "Ovarian Granulosa Cell Tumor", "Thyroid Adenoma",
               "Soft Tissue Sarcoma", "Harderian Gland Adenocarcinoma",
               "Harderian Gland Adenoma", "Osteosarcoma", "Pituitary Adenoma"), times = 3))
Tumor <- c(stack$Gamma, stack$HZE, stack$Unirradiated)


Data <- data.frame(Group, Histotype, Tumor)
Data
qplot(factor(Group), data = Data, geom = "bar", fill = factor(Tumor))

Data <- ddply(Data, .(Group),
              transform, pos = cumsum(Tumor) - (0.5 * Tumor)
)


ggplot(Data, aes(x = Group, y= Tumor, fill = Histotype),
       geom_bar(stat="identity"),
       geom_text(aes(label = Frequency, y = pos), size = 3))

p + geom_text(aes(label = Frequency), size = 3, hjust = 0.5, vjust = 3, position = "stack")




forR <- read.csv("~/Desktop/forR.csv")
forR[,"Data.Analyzed"] = forR[, "X"]


cataract <- forR[which(forR$X=="Cataract"), ]
tumor <- forR[which(forR$X=="Tumor incidence"), ]
tumor.latency <- forR[which(forR$X=="Tumor latency"), ]
survival <- forR[which(forR$X=="Survival"),]


attach(qty)
plot(Percent.Variance.Explained.2, LOD, main = "", xlab = "QTL Size", ylab = "-log10(p.value)")
plot(Size..Mb., Effect.Size, main = "", xlab = "QTL Size", ylab = "Effect Size")
plot(Percent.Variance.Explained, LOD, main = "", xlab = "Effect Size", ylab = "Percent.Variance.Explained")
abline(lm(LOD~Percent.Variance.Explained.2), col="red") # regression line (y~x) 
lines(lowess(LOD,Percent.Variance.Explained), col="blue") # lowess line (x,y)


scatterplot(Size..Mb.~Effect.Size | Data.Analyzed, data = forR, xlab = "Percent Variance Explained", 
            ylab = "QTL Confidence Interval (Mb)")

scatterplot3d(QTL.size,Percent.Variance.Explained,LOD, highlight.3d=TRUE,
              type="h",main="3D Scatterplot")


