#` This page contains code to find allelic information directly from MegaMUGA SNP
#` array.  In this way, you can determine which alleles occured in each of the 
#` individual mice, without the HMM probability information created with DOQTLs
#` genome reconstruction package. 

geno <- read.table("~/Desktop/R/QTL/extract.raw.data/GRSD/geno.txt")
total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Total-Table 1.csv")
total$corrected[which(total$Osteosarcoma > 0)]

OSA <- cbind(total$corrected[which(total$Osteosarcoma > 0)], 
      total$group[which(total$Osteosarcoma > 0)])
OSAsnps <- as.data.frame(geno[match(total$corrected[which(total$Osteosarcoma > 0)], rownames(geno)), ])
GRSDsnps <- cbind(geno[match(total[, 2], rownames(geno)), ])
GRSD.Rb1 <- cbind(GRSDsnps$UNC24228305)
OSA.Rb1 <- cbind(GRSDsnps$UNC24228305)

OSAsnps$UNC24228305 == "C"
OSAsnps$UNC24238214

rownames(Rb1)

OSA <- cbind(geno[total$corrected[which(total$Osteosarcoma == 1)], ])
NonOSA<- cbind(geno[total$corrected[which(total$Osteosarcoma == 0)], ])

snp.names <- colnames(geno)


OSA <- Total$corrected[which(Total$Osteosarcoma > 0)]
OSA <- Total$group[which(Total$Osteosarcoma > 0)]
OSA
library(HZE)
