library(knitr)
library(dplyr)
library(ggforce)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(GeomxTools)
library(readxl)
library(enrichplot)
library(data.table)
library(fgsea)
library(ggplot2)
library(ggrepel) 
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationHub)
library(GOSemSim)
library(clusterProfiler)
library(GOSemSim)
library(ggwordcloud)
library(ggplot2)
library(cowplot)
library(ReactomePA)
library(DOSE)
library(msigdbr)
library(knitr)
library(dplyr)
library(ggforce)
library(GeoMxWorkflows)
library(NanoStringNCTools)
library(GeomxTools)
library(readxl)
library(topGO)
library(scales) # for percent
library(reshape2)
library(cowplot) 
library(umap)
library(Rtsne)

load("C:/Users/edmondsonef/Desktop/DSP GeoMX/data/WTA_04122022/RData/KPC_geoMX_exp1.RData")
datadir <-"C:/Users/edmondsonef/Desktop/R-plots/"
setwd(datadir)

#test <- "trans1"
#test <- "trans2"
#test <- "trans3"
#test <- "trans4"
#test <- "trans5"
test <- "trans6"

# convert test variables to factors
pData(target_myData)$testRegion <- 
  factor(pData(target_myData)$trans6)#, c("PDAC", "Lung_met"))                           
pData(target_myData)[["slide"]] <-                                            ### Control for 
  factor(pData(target_myData)[["MHL Number"]])
assayDataElement(object = target_myData, elt = "log_q") <-
  assayDataApply(target_myData, 2, FUN = log, base = 2, elt = "q_norm")

# run LMM:
# formula follows conventions defined by the lme4 package
results <- c()
for(status in c("Full ROI")) {
  ind <- pData(target_myData)$segment == status
  mixedOutmc <-
    mixedModelDE(target_myData[, ind], elt = "log_q",
                 #modelFormula = ~ testRegion + (1 + testRegion | slide),        
                 modelFormula = ~ testRegion + (1 | slide),
                 groupVar = "testRegion",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- status
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results <- rbind(results, r_test)
}



results$Color <- "NS or FC < 0.5"
results$Color[results$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results$Color[results$FDR < 0.05] <- "FDR < 0.05"
results$Color[results$FDR < 0.001] <- "FDR < 0.001"
results$Color[abs(results$Estimate) < 0.5] <- "NS or FC < 0.5"
results$Color <- factor(results$Color, levels = c("NS or FC < 0.5", "P < 0.05", "FDR < 0.05", "FDR < 0.001"))
dplyr::count(results, FDR < 0.05)
dplyr::count(results, FDR < 0.001)

top <- dplyr::filter(results, results$FDR < 0.05)
count = count(top)
print(paste(test, ":",count, "genes with FDR < 0.05."))

head(results)
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()
for(cond in c("Full ROI")) {
  ind <- results$Subset == cond
  top_g <- c(top_g,
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = TRUE)[1:20]],
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = FALSE)[1:20]])
}
top_g <- unique(top_g)
top_g
head(top_g)



features <- c("Cybrd1","Nr1d1","Bsg","Tmprss4","Tm9sf3","Mmp23","Rhof","Sftpd", "Aqp5","Ccna1",
              "Muc3","Muc5ac","Muc3a","Kif12","Calml4","Dbp", "Mrtfb", "Rplp0","Dnajc10","Rps12",
              "Pdzd8", "Mtch2", "Msln", "Prom1", "Vars2","Porcn","Rpl6","Ybx1","Wfdc2","Tpi1",
              "Golim4","Otop3","F3", "Id2","Adamtsl5","Bag1","Rnf186","Glis2","Slc35f5","Tspan12",
              "Slc9a4", "Ephb2", "Tmem45b","Tmprss2","Pdxdc1","Lgals2", "Esrp1", "Tmem54", "Ptprf", "Ccnd2",
              "Ern2","Sult1c2", "Gltp","Spock3","Sgms2","Rasgrf1","St8sia3",
              "Rap1gap","Rbms3","Ccdc92","Ncald","Ppp1r1b","Gabbr2","Nt5c2","Cdkn2a","Atrnl1","Camk2n1",
              "Setbp1","Dennd4c","Hs3st1","Shf","Nfib", "Tuba1b", "Net1", "Ncald","Spock3",
              "Rock2", "Sem1", "Ctnnd1","Adgre5", "Dennd4c",
              "Smad4", "Flna", "Cntn1", "Cntn6","Sgms2","Nrxn1","Nrxn2","Nrxn3","Lamb2","Rasgrf1",
              "Sema3d", "Sema4b","Sema4g","Sema5a","St8sia3",
              "Lama5", "Rtn4", "Picalm","Efnb2", "Rbms3", "Rock2","Ephb2","Efnb2", "Adam10", "Mmp2", "Mmp9", 
              "Msln","Prom1","Rac1", "St8sia3", "Camk2n1", "Cdc42", "Spock3", "Rasgrf1",
              "Lama5", "Itgb1", "Ezr","S100a6", "Gsto1", "Gkn1", 
              "Lypd8l", "Anxa2", "Cdh1", "Prom1", "Myrf", "Flna", "Slc12a2", "Actn1", "Fn1", "Hnf1b",
              "Vasp","Vdac2", "Syncrip", "Rpl5", "Pard3","Dync1i2", "Calm1", "Calm2", "Calm3", "Itgb1","Kras","Trp53","Net1","Nt5c2","Ezr",
              "Clu","S100a6", "Anxa2", "Myrf", "Sema4b","Sema4g","Efnb2", "Flna", "Slc12a2", "Actn1", "Actb","Tuba1b",
              "Vasp", "Syncrip", "Pard3","Rock2","Rac1", "Rhoa", "Cdc42", "Dync1i2", "Calm1", "Calm2", "Calm3","Lama5", "Itgb1")





# #reverse log fold change to fit with label
# results$Estimate1 <- results$Estimate*(-1)
# results$Contrast[1]
# Graph results
volc_plot <- ggplot(results,                                                             ###CHANGE
       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "ALL <- log2(FC) -> Metastasis",                                       ###CHANGE
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue", `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",`NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results, Gene %in% features & `Pr(>|t|)` < 0.05),
  #geom_text_repel(data = subset(results, Gene %in% top_g & FDR < 0.01),
                  size = 6, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") 
volc_plot


ggplot(pData(target_myData), aes(x = class, y = assayDataElement(target_myData["Mmp23", ], elt = "q_norm"))) +
  geom_violin() +
  geom_jitter(width = .2) +
  labs(y = "Mmp23 Expression") +
  scale_y_continuous(trans = "log2") +
  #facet_wrap(~class) +
  theme_bw()


## ----heatmap, eval = TRUE, fig.width = 8, fig.height = 6.5, fig.wide = TRUE----
# select top significant genes based on significance, plot with pheatmap
GOI <- unique(subset(results, `FDR` < 0.01)$Gene)
pheatmap(log2(assayDataElement(target_myData[GOI, ], elt = "q_norm")),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         #clustering_method = "average",
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         cutree_cols = 3, 
         cutree_rows = 2,
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(target_myData)[, c("comps", "tissue")])

## ----maPlot, fig.width = 8, fig.height = 12, fig.wide = TRUE, warning = FALSE, message = FALSE----
results$MeanExp <-
  rowMeans(assayDataElement(target_myData,
                            elt = "q_norm"))

top_g2 <- results$Gene[results$Gene %in% top_g &
                         results$FDR < 0.01 &
                      abs(results$Estimate) > .5 &
                        results$MeanExp > quantile(results$MeanExp, 0.9)]

ggplot(subset(results, !Gene %in% neg_probes),
       aes(x = MeanExp, y = Estimate,
           size = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_hline(yintercept = c(0.5, -0.5), lty = "dashed") +
  scale_x_continuous(trans = "log2") +
  geom_point(alpha = 0.5) + 
  labs(y = "Enriched in XXX <- log2(FC) -> Enriched in XXX",
       x = "Mean Expression",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray")) +
  geom_text_repel(data = subset(results, Gene %in% features),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2) +
  theme_bw(base_size = 16) +
  facet_wrap(~Subset, nrow = 2, ncol = 1)






