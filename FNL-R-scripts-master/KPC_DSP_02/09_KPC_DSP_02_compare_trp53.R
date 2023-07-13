library(knitr)
library(dplyr)
library(ggforce)

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
load("F:/GeoMX KPC/WTA_11232022/processed_data/KPC_geoMX_exp2.RData")

#test <- "Step1_KPC"
#test <- "Step2_KPC"
#test <- "Step3_KPC_allmet"
#test <- "Step3_KPC_lung"
#test <- "Step3_KPC_liver"

#test <- "Step3_ortho"

#test <- "Step1_R172H"
#test <- "Step2_R172H"
#test <- "Step3_R172H"

#test <- "Step1_R270H"
#test <- "Step2_R270H"
#test <- "Step3_R270H"

#test <- "p53_panin"
test <- "p53_PDAC"
#test <- "p53_liver"

# convert test variables to factors
pData(target_myData)$testRegion <- 
  factor(pData(target_myData)$p53_PDAC, c("R172H_PDAC", "R270H_PDAC"))                           
pData(target_myData)[["slide"]] <-                                            ### Control for 
  factor(pData(target_myData)[["MHL"]])
assayDataElement(object = target_myData, elt = "log_q") <-
  assayDataApply(target_myData, 2, FUN = log, base = 2, elt = "q_norm")

# run LMM:
# formula follows conventions defined by the lme4 package
results <- c()
for(status in c("Full ROI")) {
  ind <- pData(target_myData)$segment == status
  mixedOutmc <-  mixedModelDE(target_myData[, ind], elt = "log_q",
                 #modelFormula = ~ testRegion + (1 + testRegion | slide),        
                 modelFormula = ~ testRegion + (1 | slide),
                 groupVar = "testRegion",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  r_test$Gene <-  unlist(lapply(colnames(mixedOutmc),
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
dplyr::count(results, `Pr(>|t|)` < 0.05)

# top <- dplyr::filter(results, `Pr(>|t|)` < 0.05)
# # top <- dplyr::filter(results, results$FDR < 0.05)
# write.csv(top, "F:/GeoMX KPC/WTA_11232022/processed_data/STRAIN_p05_KPC_PDAC.csv")

head(results)
results$invert_P <- (-log10(results$`Pr(>|t|)`)) * sign(results$Estimate)
top_g <- c()
for(cond in c("Full ROI")) {
  ind <- results$Subset == cond
  top_g <- c(top_g,
             results[ind, 'Gene'][order(results[ind, 'invert_P'], decreasing = TRUE)[1:50]],
             results[ind, 'Gene'][order(results[ind, 'invert_P'], decreasing = FALSE)[1:50]])
}
top_g <- unique(top_g)
top_g


#reverse log fold change to fit with label
results$Estimate1 <- results$Estimate*(-1)
# Graph results
volc_plot <- ggplot(results,                                                             ##
                    aes(x = Estimate1, y = -log10(`Pr(>|t|)`),
                        color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = " <- log2(FC) -> ", y = "Significance, -log10(P)", color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue", `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",`NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results, Gene %in% top_g),# & FDR < 0.01),
                  #geom_text_repel(data = subset(results, Gene %in% features & `Pr(>|t|)` < 0.05),
                  size = 6, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2, max.overlaps = 50) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") 
volc_plot


###enrichGO PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###
###enrichGO PATHWAY###

universe <- read.csv("F:/GeoMX KPC/WTA_11232022/processed_data/universePDAC.csv")
universe <- dplyr::select(universe, SYMBOL,ENTREZID)
head(universe)

R172H <- read_excel("F:/GeoMX KPC/WTA_11232022/processed_data/STRAIN_PDAC_venn_diagram.xlsx", 
                    sheet = "R172H Step2")
R270H <- read_excel("F:/GeoMX KPC/WTA_11232022/processed_data/STRAIN_PDAC_venn_diagram.xlsx", 
                    sheet = "R270H Step2")
Null <- read_excel("F:/GeoMX KPC/WTA_11232022/processed_data/STRAIN_PDAC_venn_diagram.xlsx", 
                   sheet = "Null Step2")

#### FORMAT DATA
eg <- bitr(R172H$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"),OrgDb="org.Mm.eg.db")
R172H <- dplyr::left_join(R172H, eg, by = "SYMBOL")
rm(eg)

eg <- bitr(Null$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"),OrgDb="org.Mm.eg.db")
Null <- dplyr::left_join(Null, eg, by = "SYMBOL")
rm(eg)

eg <- bitr(R270H$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"),OrgDb="org.Mm.eg.db")
R270H <- dplyr::left_join(R270H, eg, by = "SYMBOL")
rm(eg)


library(VennDiagram)
library(gridExtra)
library(readxl)
library(ggpubr)
library(Rmisc)
library(tidyverse)
library(plyr)
library(GGally)
library(ggplot2)
library(tidyverse)
library(gapminder)
library(dplyr)

Null_up <- dplyr::filter(Null, Null$Direction == "up")
R172H_up <- dplyr::filter(R172H, R172H$Direction == "up")
R270H_up <- dplyr::filter(R270H, R270H$Direction == "up")

gene_list <- list(Null = Null_up$SYMBOL, 
                  R172H = R172H_up$SYMBOL,
                  R270H = R270H_up$SYMBOL)
VennDiagram <- venn.diagram(x = gene_list, 
                            fill = c("blue", "red", "green"),
                            cat.col = c("blue", "red", "green"),
                            cex = 2,lty = "blank",
                            cat.cex = 2,
                            filename = NULL)
cowplot::plot_grid(VennDiagram)



list <- get.venn.partitions(gene_list) %>% dplyr::as_tibble()
#write.csv(list$..values..$`7`, "C:/Users/edmondsonef/Desktop/MetGenesEnrichedOverlap.csv")

list_1 <- as.data.frame(list$..values..$`2`)
colnames(list_1) <- "SYMBOL"
list_6 <- as.data.frame(list$..values..$`6`)
colnames(list_6) <- "SYMBOL"
list_4 <- as.data.frame(list$..values..$`4`)
colnames(list_4) <- "SYMBOL"
resultsGO <- dplyr::bind_rows(list_1, list_4, list_6)

colnames(resultsGO) <- "SYMBOL"
eg <- bitr(resultsGO$SYMBOL, fromType="SYMBOL", toType=c("ENTREZID"),OrgDb="org.Mm.eg.db")
resultsGO <- dplyr::left_join(resultsGO, eg, by = "SYMBOL")
rm(eg)




resultsGO <- R172H_up
summary(resultsGO)

ego <- enrichGO(gene          = resultsGO$ENTREZID,
                keyType       = "ENTREZID",
                universe      = as.character(universe$ENTREZID), 
                OrgDb         = org.Mm.eg.db,
                ont           = "BP", #"BP", "MF", "CC"
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
head(ego)

goplot(ego)
clusterprofiler::dotplot(ego)
upsetplot(ego)
#plotGOgraph(ego, useFullNames = T, useInfo = "names")


wcdf<-read.table(text=ego$GeneRatio, sep = "/")[1]
wcdf$term<-ego[,2]
wcdf$p.adjust<-ego$p.adjust
wcdf <- dplyr::filter(wcdf, V1 > 5)
wcdf <- dplyr::top_n(wcdf, 100, V1)
wcdf

selected_pathways <- c("synapse organization",
                       "axon guidance",
                       "regulation of synapse assembly",
                       "regulation of trans-synaptic signaling",
                       "positive regulation of synapse assembly",
                       "regulation of membrane potential",
                       "transmission of nerve impulse",
                       "axon development",
                       "synapse assembly",
                       "axonogenesis",
                       "chemical synaptic transmission, postsynaptic",
                       "synaptogenesis",
                       "gliogenesis",
                       "axonogenesis", 
                       "cell-substrate adhesion",
                       "oligodendrocyte development", 
                       "neurogenesis",
                       "cell-substrate adhesion",
                       "regulation of actin cytoskeleton organization",
                       "regulation of neurotransmitter transport")
dotplot(ego, showCategory = selected_pathways, font.size=10)


edox <- setReadable(ego, 'org.Mm.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=geneList)



#Create geneList
head(results)
results <- distinct(results, SYMBOL, .keep_all = T)
geneList = results$Estimate1
names(geneList) = as.character(results$ENTREZID)
geneList = sort(geneList, decreasing = T)

###gseGO PATHWAY###
ego2 <- gseGO(geneList      = geneList, 
              OrgDb        = org.Mm.eg.db,
              ont          = "BP", #"BP", "MF", and "CC"
              minGSSize    = 10,
              maxGSSize    = 1000,
              pvalueCutoff = 0.05,
              verbose      = FALSE)
goplot(ego2)
dotplot(ego2)
upsetplot(ego2, 10)




wcdf<-read.table(text=ego2$ID, sep = "/")[1]
wcdf$term<-ego2[,2]
wcdf$p.adjust<-ego2$p.adjust
wcdf <- dplyr::filter(wcdf, V1 > 20)
wcdf <- dplyr::top_n(wcdf, 100, V1)
wcdf

selected_pathways <- c("synapse organization",
                       "axon guidance",
                       "axon development",
                       "axonogenesis",
                       "synaptogenesis",
                       "gliogenesis",
                       "axonogenesis", 
                       "cell-substrate adhesion",
                       "oligodendrocyte development", 
                       "neurogenesis",
                       "cell-substrate adhesion",
                       "regulation of actin filament organization",
                       "anterograde trans-synaptic signaling",
                       "trans-synaptic signaling",
                       "synaptic signaling",
                       "regulation of dendrite extension",
                       "dendrite extension",
                       "regulation of trans-synaptic signaling")

dotplot(ego2, showCategory = selected_pathways, font.size=10)









###Compare cluster
###Compare cluster
###Compare cluster
###Compare cluster
###Compare cluster



#gcSample =  list of different samples
R172H_up[1] <- "R172H"
R270H_up[1] <- "R270H"
Null_up[1] <- "Null"
resultsGO <- dplyr::bind_rows(Null_up, R172H_up, R270H_up)
resultsCC <- dplyr::select(resultsGO, Estimate, x, ENTREZID)
unique(resultsCC$x)


#Create a list containing gene IDs for each category
ids <- unique(resultsCC$x)
mt_list<-list()
for(i in 1:length(ids)){
  id <- ids[i]
  df <- dplyr::filter(resultsCC, resultsCC$x == id)
  mt_list[[i]]<-  as.character(df$ENTREZID)
}
###CHECK THAT ROW NAMES ARE ACCURATE
names(mt_list) <- c(ids) 
str(mt_list)

ck <- compareCluster(geneCluster = mt_list, 
                     #fun = "enrichKEGG", organism = "mmu")                     #"groupGO", "enrichGO", "enrichKEGG", "enrichDO" or "enrichPathway"
                     fun = "enrichGO", ont = "BP",
                     OrgDb = org.Mm.eg.db)
ck <- setReadable(ck, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
ck <- pairwise_termsim(ck)
ck2 <- clusterProfiler::simplify(ck, cutoff=0.7, by="p.adjust", select_fun=min)

enrichplot::dotplot(ck2)
head(ck2)
goplot(ck2)
cnetplot(ck2)

wcdf<-read.table(text=ck@compareClusterResult$ID, sep = "/")[1]
wcdf$term<-ck@compareClusterResult$Description
wcdf$p.adjust<-ck@compareClusterResult$p.adjust
wcdf <- dplyr::filter(wcdf, V1 > 20)
wcdf <- dplyr::top_n(wcdf, 100, V1)
wcdf








R172H[1] <- "R172H"
R270H[1] <- "R270H"
Null[1] <- "Null"
resultsGO <- dplyr::bind_rows(Null, R172H, R270H)



resultsCC <- dplyr::select(resultsGO, Estimate, Direction, x, ENTREZID)
unique(resultsCC$x)



formula_res <- compareCluster(ENTREZID~Direction+x, 
                              data=resultsCC, 
                              fun = "enrichGO", ont = "BP",
                              OrgDb = org.Mm.eg.db)

head(formula_res)
formula_res <- clusterProfiler::simplify(formula_res, cutoff=0.7, by="p.adjust", select_fun=min)

enrichplot::dotplot(formula_res, x="Direction",
                    font.size=10, by = "Count") + facet_grid(~x)



wcdf<-read.table(text=formula_res@compareClusterResult$ID, sep = "/")[1]
wcdf$term<-formula_res@compareClusterResult$Description
wcdf$p.adjust<-formula_res@compareClusterResult$p.adjust
wcdf <- dplyr::filter(wcdf, V1 > 20)
wcdf <- dplyr::top_n(wcdf, 100, V1)
wcdf


selected_pathways <- c("angiogenesis",
                       "response to hypoxia",
                       "negative regulation of intrinsic apoptotic signaling pathway",
                       "regulation of apoptotic signaling pathway",
                       "epithelial cell apoptotic process",
                       "regulation of response to DNA damage stimulus",
                       "regulation of apoptotic signaling pathway")

enrichplot::dotplot(formula_res, x="Direction",
                    showCategory = selected_pathways, 
                    font.size=10) + facet_grid(~x)

cnetplot(formula_res, showCategory = selected_pathways)



























