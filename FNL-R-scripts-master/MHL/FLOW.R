

library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(flowAI)
library(gridExtra)
library(tidyverse)
library(flowStats)
library(flowWorkspace)
library(CytoML)
library(Rtsne)
library(FlowSOM)


myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/", pattern = ".FCS", ignore.case = TRUE)
wsp_file <- "C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/15726 10Mar2022 Simone.wsp"
wsp <- open_flowjo_xml("C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/15726 10Mar2022 Simone.wsp")
fcs_file <- "C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/Samples_Tube_018 Animal 120 BMC_018.fcs"


myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/", pattern = ".FCS", ignore.case = TRUE)
#fs <- read.flowSet(myfiles, path="C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022", truncate_max_range = FALSE)
#fs_comp <-compensate(fs, spillover(fs[[1]])$SPILL)

#tail(fj_ws_get_sample_groups(wsp))
#fj_ws_get_samples(wsp, group_id = 1)




#Removing stuff
#Removing stuff
#Removing stuff
gs <- flowjo_to_gatingset(wsp, name = 1, path ="C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/")
plot(gs)
autoplot(gs[[1]])

gs_get_pop_paths(gs)
recompute(gs)

autoplot(gs, "hCD45+")
ggcyto(gs, aes(x = `FSC-A`)) + geom_density() + geom_gate("hCD45+")

gh_pop_get_stats(gs[[5]], "hCD45+")

fs <- gh_pop_get_data(gs[[5]], "hCD45+")
fs
nrow(fs)

colnames(fs)
colnames(fs)[colnames(fs)=="Comp-BB515-A"] <- "CD8"
colnames(fs)[colnames(fs)=="Comp-BB700-P-A"] <- "CD4"
colnames(fs)[colnames(fs)=="Comp-APC-A"] <- "CD11b"
colnames(fs)[colnames(fs)=="Comp-APC-Cy7-A"] <- "CD19"
colnames(fs)[colnames(fs)=="Comp-BV421-A"] <- "CD3"
colnames(fs)[colnames(fs)=="Comp-BV786-A"] <- "CD33"
colnames(fs)[colnames(fs)=="Comp-BUV395-A"] <- "mCD45"
colnames(fs)[colnames(fs)=="Comp-BUV805-A"] <- "huCD45"
colnames(fs)[colnames(fs)=="Comp-PE-A"] <- "CD56"
colnames(fs)[colnames(fs)=="Comp-PE-CF594-A"] <- "CD66b"
colnames(fs)[colnames(fs)=="Comp-PE-Cy7-A"] <- "CD25"
colnames(fs)

fs <- flowCore::compensate(fs, flowCore::keyword(fs)[["SPILL"]])
fs <- flowCore::transform(fs,
                          flowCore::transformList(colnames(flowCore::keyword(fs)[["SPILL"]]),
                                                  flowCore::logicleTransform()))
fsom <- FlowSOM(fs, 
                #compensate = TRUE, 
                transform = TRUE,
                toTransform = c(7:12, 15, 16, 17), 
                colsToUse = c(7:12, 15, 16, 17),
                nClus = 15,
                seed = 1)
PlotStars(fsom, view = "MST", markers = fsom$map$colsUsed)
PlotStars(fsom, equalNodeSize = TRUE)

PlotDimRed(fsom, cTotal = 10000, colsToUse = fsom$map$colsUsed, 
           colorBy = "metaclusters", label = FALSE)


PlotPies(fsom, cellTypes = gatingResult$manual)






cell_types <- c("/scatter/sing/hCD45+",
                "Q6: CD3+ , CD4 [PCP55]+",
                "Q10: CD3+ , CD8 [FITC]+",
                "Q13: CD3- , CD19 [AFire750]+",
                "Q17: CD3- , CD56+",
                "Q18: CD3+ , CD56+",
                "Q29: CD66b [PEDazz]- , CD11b [AF647]+",
                "Q30: CD66b [PEDazz]+ , CD11b [AF647]+",
                "Q31: CD66b [PEDazz]+ , CD11b [AF647]-",
                "Q33: CD33- , CD11b [AF647]+",
                "Q38: CD25+ , CD3+",
                "Q35: CD33+ , CD11b [AF647]",
                "Q39: CD25+ , CD3-")

#cell_types <- c("/scatter/sing/hCD45+")

# Parse the FlowJo workspace   


gatingResult <- GetFlowJoLabels(fs, wsp_file,
                                cell_types = cell_types,
                                getData = TRUE)
# Check the number of cells assigned to each gate
colSums(gatingResult$matrix)
colnames(gatingResult$matrix)



gatingResult <- GetFlowJoLabels(fs, wsp_file,
                                cell_types = cell_types,
                                getData = TRUE)


# Build a FlowSOM tree
fsom <- FlowSOM(gatingResult$flowFrame, 
                #compensate = TRUE, 
                transform = TRUE,
                toTransform = c(7:17), 
                colsToUse = c(7:17),
                nClus = 20,
                seed = 1)
PlotStars(fsom)
PlotFlowSOM(fsom, equalNodeSize = F)

PlotPies(fsom, cellTypes = cell_types,
         backgroundValues = fsom$metaclustering)

PlotManualBars(fsom, manualVector = gatingResult$manual,
               manualOrder = c(cellTypes = cell_types))





