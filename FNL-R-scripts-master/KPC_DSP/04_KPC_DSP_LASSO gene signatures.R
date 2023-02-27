# Steps:
# 1. DE genes using DESEQ2
#     Up-regulated genes -> DAVID profling
# 2. K-means unsupervised clustering of DEGs to define molecular subtypes
# 3. Machine learning LASSO regression to get generate predictive scoring model (x number of genes, signature)
# 4. Characterize clustered with single sample get set enrichment analysis (ssGSEA)
# 
# https://jitc.bmj.com/content/10/Suppl_2/A140
# 27
# Spatial-specific gene signatures outperform bulk-mRNA signatures to define resistance to immunotherapy in melanoma patients 
#
# Abstract: Background Gene signatures have been shown to predict the response/resistance to immunotherapies but with only modest accuracy.1-3 However, gene signatures are devoid of spatial information, and the inability to distinguish tumor genes from TME genes is likely to decrease the prediction accuracy. Here we collect spatially defined genes and determine if the addition of spatial information can improve prediction.
# 
# Methods: The Digital Spatial Profiling (DSP) of GeoMX Whole Transcriptome Atlas (WTA), a new technology for transcriptome-wide spatial profiling of tissues, is used to generate a transcriptomic map of 55-immunotherapy-treated melanoma samples. DSP-WTA approach enables in situ hybridization against 18,190 genes in several areas of interest (i.e., CD68+ macrophages, CD45+ lymphocytes and S100+ tumor cells) at high throughput using a sequencing readout. We developed a computational pipeline to discover cell-type (compartment) specific signature models. To estimate the classification accuracy, the models were built using a split sample approach with 100 different Lasso logistic regression binomial models. Response Evaluation Criteria in Solid Tumors (RECIST) 1.1 was used to identify objective response. We then developed a compartment signature matrix to deconvolve bulk-RNA-seq into compartmentalized-RNA expression using CIBERSORTx (in silico tissue dissection)4 to simulate compartment-derived gene expression data. We then compare the CIBERSORT results to actual spatially collected gene signatures.
# 
# Results: We achieved AUC (area under the curve) > 0.9 for all compartment-specific signatures. The AUCs for each signature are 0.94 (CD45), 0.97 (CD68), 0.93 (S100B), 0.98 (pseudo-stroma) and 0.92 (pseudo-bulk). Cross-testing in different compartments (i.e., CD45 signature in CD68, S100B, pseudo-bulk and pseudo-stroma compartments, etc.) showed poor performance indicating the signatures are compartment-specific. Then these signatures were validated in an independent immunotherapy-treated melanoma cohort (N=90)5 by deconvolving bulk-RNA-seq into compartmentalized gene expression. Our 8-gene pseudo-bulk signature validated with AUC:0.74 (CI:0.64-0.84), our 8-gene CD68 signature with AUC:0.83 (CI:0.75-0.91), and our S100B 8-gene signature with AUC:0.66 (CI:0.54-0.77), respectively. To evaluate the deconvolution platform with real compartment RNA-seq data, our pseudo-bulk data were deconvolved into pseudo-compartments using CIBERSORTx cell-type decomposition and tested with compartment-specific signatures. We noted poor performances with AUC:0.59 (CD45), AUC:0.56 (CD68) and AUC:0.65 (S100B), respectively.
# 
# Conclusions Compartmentalized CD45-, CD68- and S100B-signatures show strikingly high AUC compared to bulk mRNA signatures. These high spatially de novo signatures far outperform the signatures that can be achieved by computational deconvolution. We believe that the spatially informed signatures differentially evaluates the tumor vs the TME and may be much more accurate in the prediction of response/resistance to immunotherapy.
# 
# https://www.biorxiv.org/content/10.1101/2021.08.01.454649v2.full.pdf
# #
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