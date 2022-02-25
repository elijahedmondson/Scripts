######EFE#######
######EFE#######
######EFE#######
######EFE#######
######EFE#######
######EFE#######
######EFE#######
######EFE#######
founder <- read.delim("C:/Users/edmondsonef/Desktop/AT_HS/AT_HS.NPT_Snps/AT_SNPs/CSU_batch1/batch1_FinalReport/batch1_FinalReport.txt", skip= 9)
write.table(Founders, file = "C:/Users/edmondsonef/Desktop/AT_HS/AT_HS.NPT_Snps/AT_SNPs/CSU_batch1/batch1_FinalReport/Founders.8xMF.1xM.17.txt", quote = FALSE, row.names = FALSE)

######EFE#######
######EFE#######
######EFE#######
######EFE#######
######EFE#######
######EFE#######
######EFE#######
######EFE#######

# parse muga data and create simplified files

library(data.table)
library(qtl2convert)

# directory containing input data
input_dir <- "../PrimaryFiles/"
# output is to MM/ and GM/

for(chip in c("MM", "GM")) {
  cat("Working on", chip, "\n")
  
  cat(" --Downloading annotations\n")
  dir <- "https://raw.githubusercontent.com/kbroman/MUGAarrays/master/UWisc/"
  ann_file <- paste0(tolower(chip), "_uwisc_v1.csv")
  if(!file.exists(ann_file)) {
    download.file(paste0(dir, ann_file), ann_file)
  }
  ann <- data.table::fread(ann_file, data.table=FALSE)
  
  cat(" --Reading data\n")
  longchip <- ifelse(chip=="MM", "MegaMUGA", "GigaMUGA")
  genofile <- paste0(input_dir, chip, "_geno.Rdata")
  snpsfile <- paste0(input_dir, chip, "_snps.Rdata")
  codefile <- paste0(input_dir, chip, "_code.Rdata")
  load(genofile)
  load(snpsfile)
  load(codefile)
  
  # simplify names
  geno <- get(paste0(chip, "_geno"))
  snps <- get(paste0(chip, "_snps"))
  code <- get(paste0(chip, "_code"))
  rm(list=c(paste0(chip, c("_geno", "_snps", "_code"))))
  
  # fix some marker names
  rownames(geno) <- gsub(".", "-", rownames(geno), fixed=TRUE)
  rownames(snps) <- gsub(".", "-", rownames(snps), fixed=TRUE)
  snps$marker <- gsub(".", "-", snps$marker, fixed=TRUE)
  
  # check that the markers are the same
  stopifnot( all(sort(rownames(geno)) == sort(ann$marker)) )
  stopifnot( all(sort(rownames(snps)) == sort(ann$marker)) )
  stopifnot( all(sort(snps$marker) == sort(ann$marker)) )
  
  # pos to Mbp
  ann$Mbp_mm10 <- ann$bp_mm10/1e6
  
  # add tier (1-4) to annotations
  ann$tier <- snps[ann$marker, "tier"]
  
  # add row names
  rownames(ann) <- ann$marker
  
  # omit markers that aren't on chr 1-19,X,Y,M
  ann <- ann[!is.na(ann$chr) & ann$chr %in% c(1:19,"X","Y","M"),]
  
  # omit non-homozygous strains
  let <- LETTERS[1:8]
  code <- code[code %in% paste0(let,let)]
  
  # omit markers not on chr 1-19,X
  # and omit the non-homozygous strains
  geno <- geno[ann$marker,names(code)]
  
  # consensus genotypes
  reduced_geno <- matrix(ncol=8, nrow=nrow(geno))
  dimnames(reduced_geno) <- list(rownames(geno), let)
  for(i in seq(along=let))
    reduced_geno[,i] <- qtl2convert::find_consensus_geno(geno[,code==paste0(let[i],let[i])])
  
  # unique alleles
  unique_alleles <- qtl2convert::find_unique_geno(reduced_geno)
  
  # omit markers that are not typed or not polymorphic
  markers2keep <- rownames(unique_alleles)[!is.na(unique_alleles[,1])]
  
  reduced_geno <- reduced_geno[markers2keep,]
  snps <- snps[markers2keep,]
  ann <- ann[markers2keep,]
  unique_alleles <- unique_alleles[markers2keep,]
  colnames(unique_alleles) <- c("A", "B")
  
  if(!dir.exists(chip)) dir.create(chip)
  
  # recode genotypes
  fg_recoded <- qtl2convert::encode_geno(reduced_geno, unique_alleles)
  
  # genetic and physical maps
  snps_gmap <- ann[,c("marker", "chr", "cM_cox")]
  colnames(snps_gmap)[3] <- "pos"
  snps_pmap <- ann[,c("marker", "chr", "Mbp_mm10")]
  colnames(snps_pmap)[3] <- "pos"
  snps_info <- ann[,c("marker", "chr", "Mbp_mm10", "cM_cox", "cM_g2f1", "tier")]
  
  # write snp info for all polymorphic markers including Y and M
  qtl2convert::write2csv(snps_info, paste0(chip, "/", chip, "_info.csv"),
                         paste("SNP information for", longchip), overwrite=TRUE)
  
  # write genotype codes
  qtl2convert::write2csv(data.frame(marker=ann$marker, chr=ann$chr, unique_alleles, stringsAsFactors=FALSE),
                         paste0(chip, "/", chip, "_allelecodes.csv"),
                         comment=paste("Allele codes for", longchip),
                         overwrite=TRUE)
  
  # write founder geno and genotype codes by chromosome
  for(chr in c(1:19,"X","Y","M")) {
    cat(" -Writing chr", chr, "\n")
    # markers on this chromosome
    mar <- ann[ann$chr==chr,1]
    
    # write founder genotypes
    qtl2convert::write2csv(data.frame(marker=mar, fg_recoded[mar,], stringsAsFactors=FALSE),
                           paste0(chip, "/", chip, "_foundergeno", chr, ".csv"),
                           comment=paste("Founder genotypes (A/B/-) for", longchip, "chromosome", chr),
                           overwrite=TRUE)
    
    # write gmap file
    qtl2convert::write2csv(snps_gmap[mar,], paste0(chip, "/", chip, "_gmap", chr, ".csv"),
                           comment=paste("Genetic map for", longchip, "chromosome", chr),
                           overwrite=TRUE)
    
    # write pmap file
    qtl2convert::write2csv(snps_pmap[mar,], paste0(chip, "/", chip, "_pmap", chr, ".csv"),
                           comment=paste("Physical map for", longchip, "chromosome", chr),
                           overwrite=TRUE)
  }
  
}