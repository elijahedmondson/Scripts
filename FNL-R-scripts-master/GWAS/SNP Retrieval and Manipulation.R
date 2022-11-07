

library(SNPtools)
snp.file = "/Users/elijah/Desktop/R/QTL/WD/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"


DOQTL:::plot.scanone.assoc(qtl, bin.size = 100)

chr = 3


#Max LOD score
top <- max(-log10(qtl[[chr]]$p.value))
top

#Position of Max LOD
max.LOD.position <- qtl[[chr]]@ranges[which(-log10(qtl[[chr]]$p.value) == top)]
max.LOD.position

max.LOD.position <- qtl[[chr]]@ranges[which(-log10(qtl[[chr]]$p.value) > (top - .1))]
max.LOD.position

start = max.LOD.position@start[1]
end = max(max.LOD.position@start[])

mgi = get.mgi.features(chr = chr, start = 54980966, end = 55200311, type = "gene", source = "MGI")
print(mgi$Name)

gene.plot(mgi)


####  SNP Tools ####

available.strains = get.strains(file = "http://cgd.jax.org/tools/SNPtools/Build38/sanger.snps.NCBI38.txt.gz")
strains = available.strains[c(4:8,11,12,14)]
strains

input = qtl[[chr]]
newqtl = data.frame(Chromosome = chr,
                    Position = input@ranges@start,
                    lod = -log10(input$p.value))


variant.plot(var.file = "http://cgd.jax.org/tools/SNPtools/Build38/sanger.snps.NCBI38.txt.gz",
        mgi.file = "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz",
        chr = chr, start = 54680966, end = 55500311, type = "snp", pattern = c("CBA/J"),
        strains = strains, ref = "CBA/J", qtl = newqtl)

snps = get.variants(chr = chr, start = start, end = end,
                    strains =  strains,
                    type = "snp")

snp = convert.variants.to.numeric(snps)

snp.plot(variants = snp, col = c("black", "grey50", "white"), cluster = T,
         ref = "C57BL/6J", highlight = "CBA/J", mgi, qtl = newqtl)


variant.plot(var.file = "http://cgd.jax.org/tools/SNPtools/Build38/sanger.snps.NCBI38.txt.gz",
             mgi.file = "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz",
             chr = chr, start = 180000000, end = 185000000, type = "snp", pattern = c("C57BL/6J"),
             strains = strains, qtl = newqtl)






layout(matrix(3:1, 3, 1))
par(mfrow = c(2,2), mar=c(1, 4, 1, 1) + 0.1)
DOQTL:::plot.scanone.assoc(HZE.days, chr=17, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Gamma.days, chr=17, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Unirradiated.days, chr=17, bin.size = 100, main = "Unirradiated", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(qtl, chr=17, bin.size = 100, main = "Total Cases", ylim=c(0,15))

par(mfrow = c(3,1), mar=c(1, 4, 1, 1) + 0.5)
DOQTL:::plot.scanone.assoc(HZE.days, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Gamma.days, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
DOQTL:::plot.scanone.assoc(Unirradiated.days, bin.size = 100, main = "Unirradiated", ylim=c(0,15))

par(mfrow = c(3,1), mar=c(1, 4, 1, 1) + 0.5)
DOQTL:::plot.scanone.assoc(HZE.days, bin.size = 100, main = "HZE Ion", ylim=c(0,15))
abline(a = 13, b = 0, col = "red")
DOQTL:::plot.scanone.assoc(Gamma.days, bin.size = 100, main = "Gamma ray", ylim=c(0,15))
abline(a = 13, b = 0, col = "red")
DOQTL:::plot.scanone.assoc(Unirradiated.days, bin.size = 100, main = "Unirradiated", ylim=c(0,15))
abline(a = 13, b = 0, col = "red")








