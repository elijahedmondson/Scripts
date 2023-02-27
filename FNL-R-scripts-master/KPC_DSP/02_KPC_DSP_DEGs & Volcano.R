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



# convert test variables to factors
pData(target_myData)$testRegion <- 
  factor(pData(target_myData)$classes, c("PDAC", "Lung_met"))                           
pData(target_myData)[["slide"]] <-                                            ### Control for 
  factor(pData(target_myData)[["Position"]])
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

head(results)
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()
for(cond in c("Full ROI")) {
  ind <- results$Subset == cond
  top_g <- c(top_g,
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = TRUE)[1:50]],
             results[ind, 'Gene'][
               order(results[ind, 'invert_P'], decreasing = FALSE)[1:50]])
}
top_g <- unique(top_g)
top_g
head(top_g)


#reverse log fold change to fit with label
results$Estimate1 <- results$Estimate*(-1)
# Graph results
ggplot(results,                                                             ###CHANGE
       aes(x = Estimate1, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = " <- log2(FC) -> ",                                       ###CHANGE
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue", `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",`NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  #geom_text_repel(data = subset(results, Gene %in% top_g & FDR < 0.01),
  geom_text_repel(data = subset(results, Gene %in% top_g),# & FDR < 0.1),
                  size = 6, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") 










ggplot(pData(target_myData),
       aes(x = classes, fill = classes,
           y = assayDataElement(target_myData["Efnb2", ],
                                elt = "q_norm"))) +
  geom_violin() +
  geom_jitter(width = .2) +
  labs(y = "") +
  scale_y_continuous(trans = "log2") +
  #facet_wrap(~Strain) +
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.title.x=element_blank(), text = element_text(size = 24))







## ----targetExprs2, fig.width = 8, fig.wide = TRUE, eval = TRUE----------------
glom <- pData(target_myData)$progression1# == "Metastasis"

# show expression of PDHA1 vs ITGB1
ggplot(pData(target_myData),
       aes(x = assayDataElement(target_myData["Rock2", ],
                                elt = "q_norm"),
           y = assayDataElement(target_myData["Efnb2", ],
                                elt = "q_norm"),
           color = classes, label=classes)) +
  geom_point(size = 3) + geom_text(hjust=1.1, vjust=0.2)+
  geom_point(size = 3) +
  theme_bw() +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  labs(x = "Rock2 Expression", y = " Expression") 
#+
#facet_wrap(~class)

## ----heatmap, eval = TRUE, fig.width = 8, fig.height = 6.5, fig.wide = TRUE----
# select top significant genes based on significance, plot with pheatmap
GOI <- unique(subset(gene, `FDR` < 0.001)$Gene)
pheatmap(log2(assayDataElement(target_myData[GOI, ], elt = "q_norm")),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         cutree_cols = 3, cutree_rows = 2,
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = pData(target_myData)[, c("progression1", "progression1")])

## ----maPlot, fig.width = 8, fig.height = 12, fig.wide = TRUE, warning = FALSE, message = FALSE----
gene$MeanExp <-
  rowMeans(assayDataElement(target_myData,
                            elt = "q_norm"))

top_g2 <- gene$Gene[gene$Gene %in% top_g &
                      gene$FDR < 0.001 &
                      abs(gene$Estimate) > .5 &
                      gene$MeanExp > quantile(gene$MeanExp, 0.9)]

ggplot(subset(gene, !Gene %in% neg_probes),
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
  geom_text_repel(data = subset(gene, Gene %in% top_g2),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2) +
  theme_bw(base_size = 16) +
  facet_wrap(~Subset, nrow = 2, ncol = 1)






