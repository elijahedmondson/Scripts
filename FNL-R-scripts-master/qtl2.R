library(qtl2)
library(devtools)
library(yaml)
library(jsonlite)
library(data.table)
library(RcppEigen)
library(RSQLite)
library(qtl)


chr <- c(1:19, "X")
write_control_file("forqtl2.json",
                   crosstype="hs",
                   description="HS/Npt AT Modifier Project",
                   founder_geno_file=paste0("GM/GM_foundergeno", chr, ".csv"),
                   founder_geno_transposed=TRUE,
                   gmap_file=paste0("GM/GM_gmap", chr, ".csv"),
                   pmap_file=paste0("GM/GM_pmap", chr, ".csv"),
                   geno_file=paste0("forqtl2_geno", chr, ".csv"),
                   geno_transposed=TRUE,
                   geno_codes=list(A=1, H=2, B=3),
                   xchr="X",
                   pheno_file="forqtl2_pheno.csv",
                   covar_file="forqtl2_covar.csv",
                   sex_covar="sex",
                   sex_codes=list(F="Female", M="Male"),
                   crossinfo_covar="ngen")