# combine info for MM and GM into a single set
# (for the case that some samples use one array and other use the other array)
#
# Of the markers on both MM and GM, there are a handful with different chromosome annotations.
# Blasting the mouse genome, these appear to map to both places equally well.
# We'll omit these markers.

# directory containing primary files
input_dir <- "PrimaryFiles"

myfread <-
    function(file)
    data.table::fread(file, data.table=FALSE)

MMinfo <- myfread("MM/MM_info.csv")
rownames(MMinfo) <- MMinfo[,1]
GMinfo <- myfread("GM/GM_info.csv")
rownames(GMinfo) <- GMinfo[,1]

cat(" -Combining info\n")

# drop markers that are in both but with inconsistent chr
MMmar <- MMinfo[,1]
GMmar <- GMinfo[,1]
in_both <- MMmar[MMmar %in% GMmar]
MMchr <- MMinfo[in_both, "chr"]
GMchr <- GMinfo[in_both, "chr"]
diff_chr <- in_both[MMchr != GMchr] # 7 markers
MMinfo <- MMinfo[!(MMinfo[,1] %in% diff_chr),]
GMinfo <- GMinfo[!(GMinfo[,1] %in% diff_chr),]

MMmar <- MMinfo[,1]
GMmar <- GMinfo[,1]
in_both <- MMmar[MMmar %in% GMmar]
MMpos <- MMinfo[in_both,"pos"]
GMpos <- GMinfo[in_both,"pos"]
names(MMpos) <- names(GMpos) <- in_both
MMcM <- MMinfo[in_both,"cM"]
GMcM <- GMinfo[in_both,"cM"]

# any missing position in one both not the other? NO
any((is.na(MMpos) & !is.na(GMpos)) | (!is.na(MMpos) & is.na(GMpos)))
any((is.na(MMcM) & !is.na(GMcM)) | (!is.na(MMcM) & is.na(GMcM)))
# cM locations are *totally* different

diff_pos <- in_both[abs(MMpos - GMpos) > 0.02]
# 47 markers that differ in position
# 20 differ by > 1 kbp
# 18 differ by >10 kbp
# I say drop those 18 and take the GM pos for the rest

# Regarding cM, I think interpolate from the GM positions

MMinfo <- MMinfo[!(MMinfo[,1] %in% diff_pos),]
GMinfo <- GMinfo[!(GMinfo[,1] %in% diff_pos),]

# combine the two
GMinfo <- cbind(GMinfo, GM=TRUE, MM=FALSE)
GMinfo[GMinfo[,1] %in% MMinfo[,1], "MM"] <- TRUE
MMinfo <- MMinfo[!(MMinfo[,1] %in% GMinfo[,1]),]
MMinfo <- cbind(MMinfo, GM=FALSE, MM=TRUE)
MMnGMinfo <- rbind(GMinfo, MMinfo)

chr <- as.numeric(factor(MMnGMinfo$chr, levels=c(1:19, "X", "Y", "M")))
pos <- MMnGMinfo$pos
pos[is.na(pos)] <- 10000 + seq(along=pos[is.na(pos)])
MMnGMinfo <- MMnGMinfo[order(chr, pos),]

for(chr in c(1:19,"X")) {
    map <- MMnGMinfo[MMnGMinfo$chr==chr & !MMnGMinfo$GM, "pos"]
    oldmap <- MMnGMinfo[MMnGMinfo$chr==chr & MMnGMinfo$GM, "pos"]
    newmap <- MMnGMinfo[MMnGMinfo$chr==chr & MMnGMinfo$GM, "cM"]
    newpos <- qtl2scan::interp_map(list("1"=map), list("1"=oldmap), list("1"=newmap))
    MMnGMinfo[MMnGMinfo$chr==chr & !MMnGMinfo$GM, "cM"] <- newpos[[1]]
}

# there were a number of markers deemed non-informative for the MegaMUGA
# but found to be informative on the GigaMUGA, and vice versa
# ...so MM and GM columns are off a bit
# correct them here
load(file.path(input_dir, "MM_snps.Rdata"))
load(file.path(input_dir, "GM_snps.Rdata"))
MM_snps[,1] <- gsub("\\.", "-", MM_snps[,1])
GM_snps[,1] <- gsub("\\.", "-", GM_snps[,1])
MMnGMinfo[!is.na(match(MMnGMinfo[,1], MM_snps[,1])),"MM"] <- TRUE
MMnGMinfo[is.na(match(MMnGMinfo[,1], MM_snps[,1])),"MM"] <- FALSE
MMnGMinfo[!is.na(match(MMnGMinfo[,1], GM_snps[,1])),"GM"] <- TRUE
MMnGMinfo[is.na(match(MMnGMinfo[,1], GM_snps[,1])),"GM"] <- FALSE

cat(" -Writing maps\n")
qtl2convert::write2csv(MMnGMinfo, "MMnGM/MMnGM_info.csv",
                       "info about combined SNPs from MegaMUGA and GigaMUGA")

for(chr in c(1:19,"X")) {

    x <- MMnGMinfo[MMnGMinfo$chr == chr, c("marker", "chr", "cM")]
    colnames(x)[3] <- "pos"
    qtl2convert::write2csv(x, paste0("MMnGM/MMnGM_", "gmap", chr, ".csv"),
                           paste("Genetic map for combined SNPs from MegaMUGA and GigaMUGA, chromosome", chr))

    x <- MMnGMinfo[MMnGMinfo$chr == chr, c("marker", "chr", "pos")]
    qtl2convert::write2csv(x, paste0("MMnGM/MMnGM_", "pmap", chr, ".csv"),
                           paste("Physical map for combined SNPs from MegaMUGA and GigaMUGA, chromosome", chr))

}

cat(" -Allele codes\n")
MMcodes <- myfread("MM/MM_allelecodes.csv")
GMcodes <- myfread("GM/GM_allelecodes.csv")
MMcodes <- MMcodes[MMcodes[,1] %in% MMnGMinfo[,1],]
GMcodes <- GMcodes[GMcodes[,1] %in% MMnGMinfo[,1],]

in_both <- MMcodes[MMcodes[,1] %in% GMcodes[,1], 1]
MMcode <- paste0(MMcodes$A, MMcodes$B)[match(in_both, MMcodes[,1])]
GMcode <- paste0(GMcodes$A, GMcodes$B)[match(in_both, GMcodes[,1])]
all(MMcode == GMcode) # yes!

MMnGMcodes <- rbind(GMcodes, MMcodes[!(MMcodes[,1] %in% GMcodes[,1]),])
MMnGMcodes <- MMnGMcodes[match(MMnGMinfo[,1], MMnGMcodes[,1]),]
qtl2convert::write2csv(MMnGMcodes, "MMnGM/MMnGM_allelecodes.csv",
                       "Allele codes for combined SNPs from MegaMUGA and GigaMUGA")

cat(" -Founder geno\n")
for(chr in c(1:19,"X")) {
    MMfg <- myfread(paste0("MM/MM_foundergeno", chr, ".csv"))
    GMfg <- myfread(paste0("GM/GM_foundergeno", chr, ".csv"))

    MMfg <- MMfg[MMfg[,1] %in% MMnGMinfo[,1],]
    GMfg <- GMfg[GMfg[,1] %in% MMnGMinfo[,1],]

    in_both <- MMfg[MMfg[,1] %in% GMfg[,1], 1]
    MMfgstr <- apply(MMfg[match(in_both, MMfg[,1]), 2:9], 1, paste, collapse="")
    GMfgstr <- apply(GMfg[match(in_both, GMfg[,1]), 2:9], 1, paste, collapse="")

    diff_fg <- in_both[MMfgstr != GMfgstr]
    x <- MMfg[match(diff_fg, MMfg[,1]),]
    y <- GMfg[match(diff_fg, GMfg[,1]),]

    # if typed in one but not the other, take the genotype
    # if typed in both but discrepant, delete genotype
    x[x != y & x=="-"] <- y[x != y & x=="-"]
    y[x != y & y=="-"] <- x[x != y & y=="-"]
    wh <- (x != y & x != "-" & y != "-")
    x[wh] <- y[wh] <- "-"

    stopifnot(all(x==y))

    # plug in differences
    MMfg[match(diff_fg, MMfg[,1]),] <- x
    GMfg[match(diff_fg, GMfg[,1]),] <- y

    MMnGMfg <- rbind(GMfg, MMfg[!(MMfg[,1] %in% GMfg[,1]), ])

    # reorder markers
    mn <- MMnGMinfo[MMnGMinfo$chr==chr,1]
    MMnGMfg <- MMnGMfg[match(mn, MMnGMfg[,1]),]

    qtl2convert::write2csv(MMnGMfg, paste0("MMnGM/MMnGM_foundergeno", chr, ".csv"),
                           "Founder genotypes (A/B/-) for combined SNPs from MegaMUGA and GigaMUGA")
}
