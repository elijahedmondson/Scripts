
Data <- ddply(Data, .(Group), 
              transform, pos = cumsum(Tumor) - (0.5 * Tumor)
)


library(plyr)
library(stringr)
library(ggplot2)
library(wesanderson)
wa.col = cbind(values=wes_palette(n=4, name="Royal1"),
               values=wes_palette(n=4, name="Zissou"),
               values=wes_palette(n=4, name="GrandBudapest2"),
               values=wes_palette(n=4, name="Royal2"))


Gamma <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Gamma-Table 1.csv")
Gamma.F <- Gamma[ which(Gamma$sex=='F'), ]
Gamma.M <- Gamma[ which(Gamma$sex=='M'), ]

Unirradiated <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/Unirradiated-Table 1.csv")
Unirradiated.F <- Unirradiated[ which(Unirradiated$sex=='F'), ]
Unirradiated.M <- Unirradiated[ which(Unirradiated$sex=='M'), ]

HZE <- read.csv("~/Desktop/R/GRSD.phenotype/CSV/HZE-Table 1.csv")
HZE.F <- HZE[ which(HZE$sex=='F'), ]
HZE.M <- HZE[ which(HZE$sex=='M'), ]


stack = data.frame(Gamma = c("Lymphoma" = sum(Gamma$Lymphoma),
                             "Myeloid Leukemia" = sum(Gamma$Myeloid.Leukemia),
                             "Pulmonary Adenocarcinoma" = sum(Gamma$Pulmonary.Adenocarcinoma),
                             "Hepatocellular Carcinoma" = sum(Gamma$Hepatocellular.Carcinoma),
                             "Hemangiosarcoma" = sum(Gamma$Hemangiosarcoma),
                             "Histiocytic Sarcoma" = sum(Gamma$Histiocytic.Sarcoma),
                             "Mammary Adenocarcinoma" = sum(Gamma$Mammary.Gland.Adenocarcinoma),
                             "Ovarian GCT" = sum(Gamma$Granulosa.Cell.Tumor),
                             "Thyroid Adenoma" = sum(Gamma$Thyroid.Tumor),
                             "Soft Tissue Sarcoma" = sum(Gamma$Soft.Tissue.Sarcomas),
                             "Harderian Gland Tumors" = sum(Gamma$Harderian.Tumor),
                             "Osteosarcoma" = sum(Gamma$Osteosarcoma),
                             "Pituitary Adenoma" = sum(Gamma$Pituitary.Adenoma)),
                   HZE = c("Lymphoma" = sum(HZE$Lymphoma),
                           "Myeloid Leukemia" = sum(HZE$Myeloid.Leukemia),
                           "Pulmonary Adenocarcinoma" = sum(HZE$Pulmonary.Adenocarcinoma),
                           "Hepatocellular Carcinoma" = sum(HZE$Hepatocellular.Carcinoma),
                           "Hemangiosarcoma" = sum(HZE$Hemangiosarcoma),
                           "Histiocytic Sarcoma" = sum(HZE$Histiocytic.Sarcoma),
                           "Mammary Adenocarcinoma" = sum(HZE$Mammary.Gland.Adenocarcinoma),
                           "Ovarian GCT" = sum(HZE$Granulosa.Cell.Tumor),
                           "Thyroid Adenoma" = sum(HZE$Thyroid.Tumor),
                           "Soft Tissue Sarcoma" = sum(HZE$Soft.Tissue.Sarcomas),
                           "Harderian Gland Tumor" = sum(HZE$Harderian.Tumor),
                           "Osteosarcoma" = sum(HZE$Osteosarcoma),
                           "Pituitary Adenoma" = sum(HZE$Pituitary.Adenoma)),
                   Unirradiated = c("Lymphoma" = sum(Unirradiated$Lymphoma),
                                    "Myeloid Leukemia" = sum(Unirradiated$Myeloid.Leukemia),
                                    "Pulmonary Adenocarcinoma" = sum(Unirradiated$Pulmonary.Adenocarcinoma),
                                    "Hepatocellular Carcinoma" = sum(Unirradiated$Hepatocellular.Carcinoma),
                                    "Hemangiosarcoma" = sum(Unirradiated$Hemangiosarcoma),
                                    "Histiocytic Sarcoma" = sum(Unirradiated$Histiocytic.Sarcoma),
                                    "Mammary Adenocarcinoma" = sum(Unirradiated$Mammary.Gland.Adenocarcinoma),
                                    "Ovarian GCT" = sum(Unirradiated$Granulosa.Cell.Tumor),
                                    "Thyroid Adenoma" = sum(Unirradiated$Thyroid.Tumor),
                                    "Soft Tissue Sarcoma" = sum(Unirradiated$Soft.Tissue.Sarcomas),
                                    "Harderian Gland Tumor" = sum(Unirradiated$Harderian.Tumor),
                                    "Osteosarcoma" = sum(Unirradiated$Osteosarcoma),
                                    "Pituitary Adenoma" = sum(Unirradiated$Pituitary.Adenoma)))

stack = data.frame(Gamma.F = c("Lymphoma" = sum(Gamma.F$Lymphoma),
                             "Myeloid Leukemia" = sum(Gamma.F$Myeloid.Leukemia),
                             "Pulmonary Adenocarcinoma" = sum(Gamma.F$Pulmonary.Adenocarcinoma),
                             "Hepatocellular Carcinoma" = sum(Gamma.F$Hepatocellular.Carcinoma),
                             "Hemangiosarcoma" = sum(Gamma.F$Hemangiosarcoma),
                             "Histiocytic Sarcoma" = sum(Gamma.F$Histiocytic.Sarcoma),
                             "Mammary Adenocarcinoma" = sum(Gamma.F$Mammary.Gland.Adenocarcinoma),
                             "Ovarian GCT" = sum(Gamma.F$Granulosa.Cell.Tumor),
                             "Thyroid Adenoma" = sum(Gamma.F$Thyroid.Tumor),
                             "Soft Tissue Sarcoma" = sum(Gamma.F$Soft.Tissue.Sarcomas),
                             "Harderian Gland Tumors" = sum(Gamma.F$Harderian.Tumor),
                             "Osteosarcoma" = sum(Gamma.F$Osteosarcoma),
                             "Pituitary Adenoma" = sum(Gamma.F$Pituitary.Adenoma)),
                   Gamma.M = c("Lymphoma" = sum(Gamma.M$Lymphoma),
                               "Myeloid Leukemia" = sum(Gamma.M$Myeloid.Leukemia),
                               "Pulmonary Adenocarcinoma" = sum(Gamma.M$Pulmonary.Adenocarcinoma),
                               "Hepatocellular Carcinoma" = sum(Gamma.M$Hepatocellular.Carcinoma),
                               "Hemangiosarcoma" = sum(Gamma.M$Hemangiosarcoma),
                               "Histiocytic Sarcoma" = sum(Gamma.M$Histiocytic.Sarcoma),
                               "Mammary Adenocarcinoma" = sum(Gamma.M$Mammary.Gland.Adenocarcinoma),
                               "Ovarian GCT" = sum(Gamma.M$Granulosa.Cell.Tumor),
                               "Thyroid Adenoma" = sum(Gamma.M$Thyroid.Tumor),
                               "Soft Tissue Sarcoma" = sum(Gamma.M$Soft.Tissue.Sarcomas),
                               "Harderian Gland Tumors" = sum(Gamma.M$Harderian.Tumor),
                               "Osteosarcoma" = sum(Gamma.M$Osteosarcoma),
                               "Pituitary Adenoma" = sum(Gamma.M$Pituitary.Adenoma)),
                   HZE.F = c("Lymphoma" = sum(HZE.F$Lymphoma),
                           "Myeloid Leukemia" = sum(HZE.F$Myeloid.Leukemia),
                           "Pulmonary Adenocarcinoma" = sum(HZE.F$Pulmonary.Adenocarcinoma),
                           "Hepatocellular Carcinoma" = sum(HZE.F$Hepatocellular.Carcinoma),
                           "Hemangiosarcoma" = sum(HZE.F$Hemangiosarcoma),
                           "Histiocytic Sarcoma" = sum(HZE.F$Histiocytic.Sarcoma),
                           "Mammary Adenocarcinoma" = sum(HZE.F$Mammary.Gland.Adenocarcinoma),
                           "Ovarian GCT" = sum(HZE.F$Granulosa.Cell.Tumor),
                           "Thyroid Adenoma" = sum(HZE.F$Thyroid.Tumor),
                           "Soft Tissue Sarcoma" = sum(HZE.F$Soft.Tissue.Sarcomas),
                           "Harderian Gland Tumor" = sum(HZE.F$Harderian.Tumor),
                           "Osteosarcoma" = sum(HZE.F$Osteosarcoma),
                           "Pituitary Adenoma" = sum(HZE.F$Pituitary.Adenoma)),
                   HZE.M = c("Lymphoma" = sum(HZE.M$Lymphoma),
                           "Myeloid Leukemia" = sum(HZE.M$Myeloid.Leukemia),
                           "Pulmonary Adenocarcinoma" = sum(HZE.M$Pulmonary.Adenocarcinoma),
                           "Hepatocellular Carcinoma" = sum(HZE.M$Hepatocellular.Carcinoma),
                           "Hemangiosarcoma" = sum(HZE.M$Hemangiosarcoma),
                           "Histiocytic Sarcoma" = sum(HZE.M$Histiocytic.Sarcoma),
                           "Mammary Adenocarcinoma" = sum(HZE.M$Mammary.Gland.Adenocarcinoma),
                           "Ovarian GCT" = sum(HZE.M$Granulosa.Cell.Tumor),
                           "Thyroid Adenoma" = sum(HZE.M$Thyroid.Tumor),
                           "Soft Tissue Sarcoma" = sum(HZE.M$Soft.Tissue.Sarcomas),
                           "Harderian Gland Tumor" = sum(HZE.M$Harderian.Tumor),
                           "Osteosarcoma" = sum(HZE.M$Osteosarcoma),
                           "Pituitary Adenoma" = sum(HZE.M$Pituitary.Adenoma)),
                   Unirradiated.F = c("Lymphoma" = sum(Unirradiated.F$Lymphoma),
                                    "Myeloid Leukemia" = sum(Unirradiated.F$Myeloid.Leukemia),
                                    "Pulmonary Adenocarcinoma" = sum(Unirradiated.F$Pulmonary.Adenocarcinoma),
                                    "Hepatocellular Carcinoma" = sum(Unirradiated.F$Hepatocellular.Carcinoma),
                                    "Hemangiosarcoma" = sum(Unirradiated.F$Hemangiosarcoma),
                                    "Histiocytic Sarcoma" = sum(Unirradiated.F$Histiocytic.Sarcoma),
                                    "Mammary Adenocarcinoma" = sum(Unirradiated.F$Mammary.Gland.Adenocarcinoma),
                                    "Ovarian GCT" = sum(Unirradiated.F$Granulosa.Cell.Tumor),
                                    "Thyroid Adenoma" = sum(Unirradiated.F$Thyroid.Tumor),
                                    "Soft Tissue Sarcoma" = sum(Unirradiated.F$Soft.Tissue.Sarcomas),
                                    "Harderian Gland Tumor" = sum(Unirradiated.F$Harderian.Tumor),
                                    "Osteosarcoma" = sum(Unirradiated.F$Osteosarcoma),
                                    "Pituitary Adenoma" = sum(Unirradiated.F$Pituitary.Adenoma)),
                   Unirradiated.M = c("Lymphoma" = sum(Unirradiated.M$Lymphoma),
                                    "Myeloid Leukemia" = sum(Unirradiated.M$Myeloid.Leukemia),
                                    "Pulmonary Adenocarcinoma" = sum(Unirradiated.M$Pulmonary.Adenocarcinoma),
                                    "Hepatocellular Carcinoma" = sum(Unirradiated.M$Hepatocellular.Carcinoma),
                                    "Hemangiosarcoma" = sum(Unirradiated.M$Hemangiosarcoma),
                                    "Histiocytic Sarcoma" = sum(Unirradiated.M$Histiocytic.Sarcoma),
                                    "Mammary Adenocarcinoma" = sum(Unirradiated.M$Mammary.Gland.Adenocarcinoma),
                                    "Ovarian GCT" = sum(Unirradiated.M$Granulosa.Cell.Tumor),
                                    "Thyroid Adenoma" = sum(Unirradiated.M$Thyroid.Tumor),
                                    "Soft Tissue Sarcoma" = sum(Unirradiated.M$Soft.Tissue.Sarcomas),
                                    "Harderian Gland Tumor" = sum(Unirradiated.M$Harderian.Tumor),
                                    "Osteosarcoma" = sum(Unirradiated.M$Osteosarcoma),
                                    "Pituitary Adenoma" = sum(Unirradiated.M$Pituitary.Adenoma)))


Group <- c(rep(c("Gamma Female", "Gamma Male", "HZE Female", "HZE Male", "Unirradiated Female", "Unirradiated Male"), each = 13))
Histotype <- c(rep(c("Lymphoma", "Myeloid Leukemia", 
                     "Pulmonary Adenocarcinoma", "Hepatocellular Carcinoma", 
                     "Hemangiosarcoma", "Histiocytic Sarcoma", 
                     "Mammary Adenocarcinoma",
                     "Ovarian GCT", "Thyroid Adenoma", 
                     "Soft Tissue Sarcoma", "Harderian Gland Tumor",
                     "Osteosarcoma", "Pituitary Adenoma"), times = 6))
Tumor <- c(stack$Gamma.F, stack$Gamma.M, stack$HZE.F, stack$HZE.M, stack$Unirradiated.F, stack$Unirradiated.M)

Data <- data.frame(Group, Histotype, Tumor)
Data
Data <- ddply(Data, .(Group), transform, pos = cumsum(Tumor) - (0.5 * Tumor))


Data$Histotype <- factor(Data$Histotype, levels = c("Lymphoma", 
                                                    "Myeloid Leukemia", "Pulmonary Adenocarcinoma", 
                                                    "Hepatocellular Carcinoma", "Hemangiosarcoma", 
                                                    "Histiocytic Sarcoma", "Mammary Adenocarcinoma", 
                                                    "Ovarian GCT", "Thyroid Adenoma", 
                                                    "Soft Tissue Sarcoma", "Harderian Gland Tumor", 
                                                    "Osteosarcoma", "Pituitary Adenoma"))
Data$Histotype <- factor(Data$Histotype, levels = rev(levels(Data$Histotype)))


#Wraps the Text for labels
Data$NewHistotype <- str_wrap(Data$Histotype, width = 15)
Data$NewGroup <- str_wrap(Data$Group, width = 5)

ggplot(Data, aes(x = NewGroup, y= Tumor)) + 
       geom_bar(aes(fill = Histotype), stat="identity", colour = "black") +
        #scale_y_continuous() +
       geom_text(aes(label = NewHistotype, y = pos, size = Tumor)) + 
        scale_radius(range = c(2,8)) +
        scale_fill_manual(values = wa.col) +
        ylab("Number of Mice") + xlab("") +
        theme_bw(base_size = 25)


                             


p + geom_text(aes(label = Frequency), size = 3, hjust = 0.5, vjust = 3, position = "stack") 

library(RColorBrewer)


my.col <- c("#ccf2ff", "#9CB071", "#9CB071",
            "#9CB071", "#9CB071", "#87AFC7",
            "#FFCBA4", "#ff9966", "#E6E600FF",
            "#2B547E", "#E8C034FF", "#404040", 
            "#e6e6e6", "#ff6666", "#617C58", 
            "#87CEEB", "#C48793", "#EDC9AF")
barplot(t(stack), col = my.col, ylab = "Number of Mice", ylim = c(0,800), xlim = c(0,12), width = 2)


legend("bottomright", 
       legend = c(colnames(stack)[18:1]), #in order from top to bottom
       fill = my.col[18:1], # 6:1 reorders so legend order matches graph
       title = "Tumor Histotype")