## Determining haplotype
load("/Users/elijah/Desktop/R/EC2 Cloud/build/markers.Rdata")
load("/Users/elijah/Desktop/R/QTL/WD/hs.colors.Rdata")

Total <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/GRSD.pheno.csv")
pheno = data.frame(row.names = Total$row.names, rownames = Total$row.names,
                   sex = as.numeric(Total$sex == "M"),
                   group = as.character(Total$groups),
                   days = as.numeric(Total$days),
                   AML = as.numeric(Total$Myeloid.Leukemia),
                   AML.ASXLdel = as.numeric(Total$Asxl1.Deletion),
                   AML.PU.1del = as.numeric(Total$Pu.1.Deletion))
AML = pheno[which(pheno$AML == 1),]
rm(pheno, Total)
#AML.ASXLdel = pheno[which(pheno$AML.ASXLdel == 1),]
#AML.PU.1del = pheno[which(pheno$AML.PU.1del == 1),]

get.allele(chr = 2, pos = 147914389, markers, pheno = AML)

allele = as.data.frame(ALLELE)

hist(allele$`allele beginning`)

alle = transform(allele, allele$`allele beginning` = colsplit(allele$`allele beginning`, split = "", names = c('a', 'b')))



get.allele = function(chr, pos, markers, pheno) {
        
        ALLELE = matrix(0, nrow = nrow(pheno), ncol = 11, dimnames = list(1:nrow(pheno), 
                        c("Mouse", "sex", "exposure", "days", "AML", "Asxl1 deletion", 
                          "PU.1 deletion", "allele beginning", "allele end", "Marker min", "Marker max")))
        
        for(i in 1:nrow(pheno)) {
                tryCatch({
                        file = paste0("~/Desktop/R/QTL/HMM/", pheno[i,1], ".genotype.probs.Rdata")
                        load(file = file)
                        CASE = markers[which(markers$Chr == chr),]
                        CASE = CASE[which(CASE$Mb_NCBI38 > (pos - 50000)),]
                        CASE = CASE[which(CASE$Mb_NCBI38 < (pos + 50000)),]
                        
                        samples = intersect(CASE$SNP_ID, rownames(prsmth))
                        stopifnot(length(samples) > 0)
                        CASE = CASE[samples,,drop = FALSE]
                        
                        prsmth = prsmth[samples,,drop = FALSE]
                        allele = colnames(prsmth)[max.col(prsmth)]
                        
                        ALLELE[i,] = c(pheno[i,1],pheno[i,2],pheno[i,3],pheno[i,4],pheno[i,5],pheno[i,6],
                                       pheno[i,7], allele[1], max(allele), (pos - 50000), (pos + 50000))
                        print(ALLELE[i,])
                }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
                
        }
        write.csv(ALLELE, file = paste0("~/Desktop/ALLELE", ".CHR = ", chr, ".POS = ", pos, ".csv"))
        print(ALLELE)
       
} # get.allele()


