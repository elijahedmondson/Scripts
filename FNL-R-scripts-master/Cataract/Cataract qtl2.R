library(qtl2)
library(qtl2convert)
library(yaml)
library(jsonlite)
library(data.table)
library(RcppEigen)
library(RSQLite)
library(qtl)
library(qtl2convert)
library(readxl)
library(dplyr)


# pheno <- read_excel("C:/Users/edmondsonef/Desktop/Cataract/CATARACT_final.xlsx", sheet ="CATARACT_fin")
# #load("C:/Users/edmondsonef/Desktop/QTL/CAT_QTLproject_2022.Rdata")
# load("C:/Users/edmondsonef/Desktop/QTL/HZEproject.RData")
# #load(url("ftp://ftp.jax.org/MUGA/MM_snps.Rdata")) 
# #load("C:/Users/edmondsonef/Desktop/QTL/Rqtl2 Files/HZE-001_qtl2.RData")
# 
# 
# pheno = data.frame(row.names = Total$row.names, sex = as.numeric(Total$sex == "F"),
#                    albino = as.numeric(Total$`coat color`=="albino"),
#                    cat = as.numeric(Total$cat_score),
#                    days = as.numeric(Total$days))
# addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))
# 
# 
# 
# 
# library(qtl2convert)
# library(qtl2)
# 
# base_dir  = 'C:/Users/edmondsonef/Desktop/QTL/Rqtl2 Files/'
# #src_file  = file.path(base_dir, 'HZE-001.Rdata')
# dest_file = file.path(base_dir, 'HZE-EFE_qtl2.Rdata')
# 
# #load(src_file)
# 
# new_probs = qtl2convert::probs_doqtl_to_qtl2(probs = probs, map = markers,
#                                              chr_column = 'Chr', pos_column = 'cM',
#                                              marker_column = 'SNP_ID')
# 
# K = calc_kinship(probs = new_probs, type = 'loco')
# 
# save(pheno, new_probs, K, markers, file = dest_file)

# map = qtl2convert::map_df_to_list(markers, chr_column = "Chr", pos_column = "Mb_NCBI38", marker_column = "SNP_ID")
# 
# total <- read_excel("C:/Users/edmondsonef/Desktop/Cataract/CATARACT_final.xlsx", sheet ="simple")
# pheno <- data.frame(row.names = total$row.names, 
#                     sex = as.numeric(total$sex == "M"), 
#                     group = as.character(total$groups),
#                     albino = as.numeric(total$`coat color` == "albino"),
#                     cat_average = as.numeric(total$`cat_average`),
#                     cat_sum = as.numeric(total$`cat_sum`),
#                     cat_difference = as.numeric(total$`cat_difference`)) 
#
#save(total, pheno, probs, K, markers, map, addcovar, file = "C:/Users/edmondsonef/Desktop/QTL/Rqtl2 Files/HZE-EFE_qtl2_cat.RData")




load("C:/Users/edmondsonef/Desktop/QTL/Rqtl2 Files/HZE-EFE_qtl2_cat.RData")

HZE <- subset(pheno, group == "HZE")
Gamma <- subset(pheno, group == "Gamma")
Un <- subset(pheno, group == "Unirradiated")
All_irr <- subset(pheno, group != "Unirradiated")


addcovar = matrix(pheno$sex, ncol = 1, dimnames = list(rownames(pheno), "sex"))
#pheno_a = data.frame(row.names = row.names(pheno), albino = as.numeric(pheno$`albino`)) 
HZE_c = data.frame(row.names = row.names(HZE), cat = as.numeric(HZE$`cat_average`))   

head(total)
head(pheno)
head(pheno_c)


###Run QTL Analysis
all_cataract <- scan1(genoprobs = probs, 
                    pheno = HZE_c, 
                    kinship = K, 
                    addcovar = addcovar, 
                    #Xcovar = NULL, 
                    #intcovar = NULL, 
                    #weights = NULL, 
                    #reml = TRUE,
                    #model = c("normal", "binary"), 
                    #hsq = NULL,
                    cores=3)


ymx <- maxlod(all_cataract) # overall maximum LOD score
find_peaks(all_cataract, map= map, threshold=4, drop=1.5)

plot(all_cataract, map, lodcolumn=1, ylim=c(0, ymx*1.02))



operm <- scan1perm(probs, 
                   pheno = HZE_c,
                   addcovar = addcovar, 
                   n_perm=1000,
                   perm_Xsp=TRUE, 
                   chr_lengths=chr_lengths(map))
summary(operm)
summary(operm, alpha=c(0.2, 0.05))




out_gwas <- scan1snps(genoprobs = probs,
                      pheno = HZE_c, 
                      kinship = K, 
                      map, 
                      #query_func=query_variants, 
                      cores=3)



