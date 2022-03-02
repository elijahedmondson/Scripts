#' # Identify the files
#' 
#' "C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15701 02Feb2022 Simone/"
#'
myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15701 02Feb2022 Simone/", pattern = ".FCS", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15701 02Feb2022 Simone/")#, truncate_max_range = FALSE)

fcs_file <- system.file("C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15716 28Feb2022/", "--------.fcs", package = "FlowSOM")
wsp_file <- system.file("C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15716 28Feb2022/", "15716 28Feb2022 Simone.wsp", package = "FlowSOM")

# Specify the cell types of interest for assigning one label per cell
cell_types <- c("B cells", "T cells", "CD4 T cells", "CD8 T cells",
                "NK cells","NK T cells")

# Parse the FlowJo workspace   
library(flowWorkspace)             
gatingResult <- GetFlowJoLabels(fcs_file, wsp_file,
                                cell_types = cell_types)

# Check the number of cells assigned to each gate
colSums(gatingResult$matrix)

# Build a FlowSOM tree
flowSOM.res <- FlowSOM(fcs_file, 
                       compensate = TRUE, 
                       transform = TRUE,
                       toTransform = 8:18, 
                       colsToUse = c(9,12,14:18),
                       nClus = 10,
                       seed = 1)

# Plot pies indicating the percentage of cell types present in the nodes
PlotPies(flowSOM.res$FlowSOM,
         gatingResult$manual,
         backgroundValues = flowSOM.res$metaclustering)