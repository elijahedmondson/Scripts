library(DOQTL)
setwd("/hpcdata/dgatti/HS/")
load(file = "HS_mapping_data.Rdata")
sdp.file = "/hpcdata/dgatti/HS_Sanger_SDPs.txt.bgz"

pheno[,2] = as.numeric(pheno[,2]) - 1

qtl = scanone.assoc(pheno = pheno, pheno.col = 5, probs = probs, K = K, 
      addcovar = pheno[,2,drop = F], markers = snps, sdp.file = sdp.file, ncl = 20)
save(qtl, file = "albino_assoc_QTL.Rdata")

png("albino_QTL.png", width = 2000, height = 1600, res = 128)
DOQTL:::plot.scanone.assoc(qtl, bin.size = 100)
dev.off()

png("albino_QTL_chr7.png", width = 2000, height = 1600, res = 128)
DOQTL:::plot.scanone.assoc(qtl, chr = 7, bin.size = 10)
dev.off()
