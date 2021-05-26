library(flowViz)
data(GvHD)
head(pData(GvHD))



library(flowStats)
browseVignettes("flowStats")

library(flowClust)
browseVignettes("flowClust")

###                            
install.packages("devtools")
devtools::install_github("r-lib/pillar")
###

library(flowCore)
library(ggcyto)

browseVignettes("flowCore")
file.name <- system.file("extdata","0877408774.B08",
                         package="flowCore")
x <- read.FCS('C:/Users/edmondsonef/Desktop/Humanized Mice/Flow Data/NSGS spleen/Samples_Tube_020 spleen 68_008.fcs', transformation=FALSE)
summary(x)
fs <- transform(x, transformList(c("FSC-A", "FSC-H"), list(log, log)))
autoplot(x, "FSC-A")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("flowCore")



# load your omic data here as mydata
library(M3C)
browseVignettes("M3C")
data(x)
tsne(x$data,colvec=c('gold'))
