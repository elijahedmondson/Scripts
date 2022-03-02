
#https://rpubs.com/artur_matysik/flow
###https://jchellmuth.com/posts/FACS-with-R/

library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(flowViz)
library(flowStats)
library(scales)
library(dplyr)
library(stringr)
library(viridis)
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(flowAI)
library(gridExtra)
library(tidyverse)
library(flowStats)

#manual
#Load data
myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15701 02Feb2022 Simone/", pattern = ".FCS", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15701 02Feb2022 Simone/")#, truncate_max_range = FALSE)
myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/samp/", pattern = ".FCS", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="C:/Users/edmondsonef/Desktop/samp/", truncate_max_range = FALSE)
fs_comp <-compensate(fs, spillover(fs[[1]])$SPILL)





pData(fs) %>% head(3)
pData(fs)$well <- gsub(".*_.*_(.*)_.*.fcs","\\1",sampleNames(fs)) # extract well from name and add new 'well' column
pData(fs) %>% head(3)
colnames(fs)

colnames(fs)[colnames(fs)=="BB515-A"] <- "CD8"
colnames(fs)[colnames(fs)=="BB700-P-A"] <- "CD4"
colnames(fs)[colnames(fs)=="APC-A"] <- "CD11b"
colnames(fs)[colnames(fs)=="APC-Cy7-A"] <- "CD19"
colnames(fs)[colnames(fs)=="BV421-A"] <- "CD3"
colnames(fs)[colnames(fs)=="BV786-A"] <- "CD33"
colnames(fs)[colnames(fs)=="BUV395-A"] <- "mCD45"
colnames(fs)[colnames(fs)=="BUV805-A"] <- "huCD45"
colnames(fs)[colnames(fs)=="PE-A"] <- "CD56"
colnames(fs)[colnames(fs)=="PE-CF594-A"] <- "CD66b"
colnames(fs)[colnames(fs)=="PE-Cy7-A"] <- "CD25"
#fs <- fsApply(fs, function(x) {transform(x, estimateLogicle(x, c("CD8", "CD4", "CD11b","CD19", "CD3", "CD33",
#                                                                       "mCD45", "huCD45", "CD56", "CD66b")))})
#fs <- fs_trans
gs <- GatingSet(fs)
g.singlets <- polygonGate(filterId = "Singlets","FSC-A"=c(2e4,26e4,26e4,2e4),"FSC-H"=c(0e4,16e4,24e4,6e4)) # define gate
gs_pop_add(gs,g.singlets) # add gate to GatingSet

ggcyto(fs[[1]], aes(x="FSC-A",y="FSC-H"),subset="root") + geom_hex(bins = 500) + 
  labs(title = "huCD45 vs FSC", y = "FSC-H", x = "FSC-A") + ggcyto_par_set(limits = "instrument") +
  geom_gate(g.singlets) +
  theme_bw() + theme(legend.position = "none", aspect.ratio = 1) + facet_wrap(~name, ncol = 4) 
gs_pop_add(gs,g.singlets) # add gate to GatingSet


recompute(gs)

ggcyto(gs[[1]],aes(x="FSC-A",y="SSC-A"),subset="root")+geom_hex(bins = 200)+ggcyto_par_set(limits = "instrument")
ggcyto(gs[[1]],aes(x="FSC-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 200)+ggcyto_par_set(limits = "instrument")

ggcyto(gs,aes(x="FSC-A",y="FSC-H"),subset="root")+geom_hex(bins = 100)+geom_gate("Singlets")+
  geom_stats(adjust = 0.8)+ggcyto_par_set(limits = "instrument")+
  facet_wrap(~well,ncol = 10)


g.live <- polygonGate(filterId = "Live","FSC-A"=c(8e4,28e4,28e4,8e4),"SSC-A"=c(3e4,4e4,28e4,28e4)) # define gate
ggcyto(gs[[1]],aes(x="FSC-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 200)+geom_gate(g.live)+ggcyto_par_set(limits = "instrument") # check gate
gs_pop_add(gs,g.live,parent="Singlets") # add gate to GatingSet
recompute(gs)


ggcyto(gs,aes(x="FSC-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 100)+geom_gate("Live")+
  geom_stats(adjust = 0.8)+ggcyto_par_set(limits = "instrument")+
  facet_wrap(~well,ncol = 10)

g.huCD45 <- rectangleGate(filterId="huCD45 positive","huCD45"=c(1000, Inf)) # set gate
ggcyto(gs[[1]],aes(x=huCD45),subset="Live")+geom_density(fill="forestgreen")+geom_gate(g.huCD45)+ggcyto_par_set(limits = "instrument")+scale_x_flowJo_biexp() # check gate
gs_pop_add(gs,g.huCD45,parent="Live")
recompute(gs)


ggcyto(gs,aes(x=huCD45),subset="Live",)+geom_density(fill="forestgreen")+geom_gate("huCD45 positive")+
  geom_stats()+
  ggcyto_par_set(limits = "instrument")+scale_x_flowJo_biexp()+
  facet_wrap(~well,ncol = 10)

gs_pop_get_count_fast(gs) %>% head
ps <- data.frame(gs_pop_get_count_fast(gs))

ps$percent_of_parent <- ps$Count/ps$ParentCount*100
psm <- merge(ps,pData(fs),by="name")




library(flowCore)
library(FlowSOM)
library(ggplot2)

# 18. Compute the FlowSOM object

SOM_x <- 10
SOM_y <- 10
n_meta <- 8
seed <- 2020
scaling <- FALSE


names(fs)
exprs(fs)
markers_of_interest <- c("SSC-A", "CD8", "CD4", "CD3", "CD33",
                         "CD11b", "CD25", "mCD45", "huCD45", "CD56", 
                         "CD19", "CD66b")
markers_of_interest <- c("SSC-A", "CD4", "CD3", "CD33",
                         "CD11b", "CD25", "mCD45", "huCD45", "CD56", 
                         "CD19", "CD66b")


fsom <- FlowSOM(input = fs,
                scale = scaling,
                colsToUse = c(4,7:16),
                seed = seed,
                nClus = n_meta,
                xdim = SOM_x, ydim = SOM_y)
saveRDS(fsom, paste(dir_RDS, "fsom.rds"))

PlotDimRed(fsom, colsToUse = fsom$map$colsUsed, colorBy = "metaclusters")


