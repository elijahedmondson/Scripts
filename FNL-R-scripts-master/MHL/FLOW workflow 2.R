library(flowCore)
library(Rtsne)
library(FlowSOMworkshop)
library(FlowSOM)

wsp_file <- "C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/15726 10Mar2022 Simone.wsp"
fcs_file <- "C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/Samples_Tube_015 Animal 137 blood_015.fcs"
tube_13 <- parse_flowjo(fcs_file, wsp_file)

fs <- flowCore::exprs(flowCore::read.FCS(fcs_file, transformation = FALSE, truncate_max_range = FALSE))

head(fs)
dim(fs)
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
head(fs)
dim(fs)

live <- gating_subset(tube_13, "Live")
live_ff <- live$flowFrame