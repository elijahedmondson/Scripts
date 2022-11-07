library(DOQTL)
library(QTLRel)
load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata"))            
                 
extract.raw.data(in.path = "/Volumes/External Terabyte/QTL/Founders", prefix = "",
                 out.path = "/Volumes/External Terabyte/QTL/extract.raw.data/Founders", 
                 array = "megamuga")

extract.raw.data(in.path = c("/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/1 - 96",
                             "/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/2 - 23 (last plate)",
                             "/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/2 - 481",
                             "/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/3 - 18 (last plate)",
                             "/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/4 - 600",
                             "/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/5 - 12 (last plate 600)",
                             "/Users/elijahedmondson/Desktop/R/QTL/FINAL GENOTYPE DATA/6 - 648"), 
                 prefix = c("", "", "", "", "", "", ""),
                 out.path = "/Users/elijahedmondson/Desktop/R/QTL/extract.raw.data/GRSD", 
                 array = "megamuga")
