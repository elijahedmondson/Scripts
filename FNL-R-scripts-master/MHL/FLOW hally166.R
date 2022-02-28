

 if (!requireNamespace("BiocManager", quietly = TRUE))
   install.packages("BiocManager")
 BiocManager::install("flowAI")

library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(flowAI)
library(gridExtra)

#manual
#Load data
myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15709 15Feb2022 Simone/", pattern = ".FCS", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15709 15Feb2022 Simone/")
fs_comp <-compensate(fs, spillover(fs[[1]])$SPILL)
fs_comp_clean <- flow_auto_qc(fs_comp)
trans <- estimateLogicle(fs_comp_clean[[1]], colnames(fs_comp_clean[,3:7]))
fs_comp_clean_trans <- transform(fs_comp_clean, trans)

#Visualise file
fs_comp_clean_trans[[1]]
autoplot(fs_comp_clean_trans[[1]])

#Basic gating
ggcyto(fs_comp_clean_trans[[1]], aes(x="FSC-A", y="SSC-A"))+geom_hex(bins=256)

#create the empty gating set
gs<-GatingSet(fs_comp_clean_trans)

#Cells - FSC SSC
rg1<-rectangleGate("FSC-A"=c(15000, Inf), filterId = "NoneDebris")
gs_pop_add(gs, rg1, parent="root")
recompute(gs)
gs_get_pop_paths(gs)
ggcyto(fs_comp_clean_trans[[1]], aes(x="FSC-A", y="SSC-A"))+geom_hex(bins=256)+geom_gate(gs_pop_get_gate(gs, "NoneDebris"))
gs_pop_get_stats(gs)

#Singlet gating
ggcyto(fs_comp_clean_trans[[1]], aes(x = "FSC-H", y = 'FSC-W'))+ geom_hex(bins = 256)
rg2 <- rectangleGate("FSC-H"=c(3.6, 4.2),"FSC-W"=c(50, 240))
gs_pop_add(gs, rg2, parent = "NoneDebris", name = "singlets")
gs_get_pop_paths(gs)
recompute(gs)
ggcyto(fs_comp_clean_trans, aes(x = "FSC-H", y = 'FSC-W'))+ geom_hex(bins = 256)+ geom_gate(gs_pop_get_gate(gs, "singlets"))

#exploring the gatingset
plot(gs)
gs_pop_get_stats(gs)
gs_pop_get_stats(gs, "NoneDebris", "percent")


#automatic
#Load data
myfiles <- list.files(path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZV_files", pattern = ".FCS", ignore.case = TRUE)
fs <- read.flowSet(myfiles[1:2], path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZV_files", alter.names=TRUE)
fs_comp <-compensate(fs,spillover(fs[[1]])$SPILL)
#fix compensation matrix
matrix<-spillover(fs[[1]])$SPILL
matrix
colnames(matrix)<-c("X.FITC.A.", "X.Pacific.Blue.A.", "X.Alexa.680.A.", "X.APC.A.", "X.PE.Cy7.A.", "X.PE.Cy55.A.", "X.PE.Tx.RD.A.", "X.PE.Green.laser.A.")
fs_comp <-compensate(fs,matrix)
#continue
fs_comp_clean <- flow_auto_qc(fs_comp)
trans <- estimateLogicle(fs_comp_clean[[1]], colnames(fs_comp_clean[,c(4,6:12)]))
fs_comp_clean_trans <- transform(fs_comp_clean, trans)
autoplot(fs_comp_clean_trans[[1]])

#create the empty gating set
auto_gs<-GatingSet(fs_comp_clean_trans)

#cell gate
fs_data<- gs_pop_get_data(auto_gs)
noneDebris_gate<- fsApply(fs_data, function(fr) openCyto:::.flowClust.2d(fr, channels= c("FSC.A","SSC.A")))
gs_pop_add(auto_gs, noneDebris_gate, parent = "root", name="noneDebris_gate")
recompute(auto_gs)
autoplot(auto_gs, x="FSC.A", y="SSC.A", "noneDebris_gate", bins=256)

#Singlet gate
fs_data <- gs_pop_get_data(auto_gs, "noneDebris_gate") #get parent data
singlet_gate <- fsApply(fs_data, function(fr) openCyto:::.singletGate(fr, channels =c("FSC.A", "FSC.H")))
gs_pop_add(auto_gs, singlet_gate, parent = "noneDebris_gate", name = "singlets")
recompute(auto_gs)
autoplot(auto_gs, x = 'FSC.A', y = 'FSC.H', "singlets", bins = 256)

#Quad gate
fs_data <- gs_pop_get_data(auto_gs, "singlets") #get parent data
BGquad_gate <- fsApply(fs_data, function(fr) openCyto:::.quadGate.seq(fr, gFunc="mindensity", min=c(3,3), channels =c('X.FITC.A.', 'X.PE.Tx.RD.A.')))
gs_pop_add(auto_gs, BGquad_gate, parent = "singlets", names = c("1", "2", "3", "4"))
recompute(auto_gs)
gs_get_pop_paths(auto_gs[[1]])
plot(auto_gs)
autoplot(auto_gs, x = 'X.FITC.A.', y = 'X.PE.Tx.RD.A.', gs_get_pop_paths(auto_gs)[4:7], bins = 256)

#fix plot
p<-ggcyto(auto_gs[1:2],aes(x = 'X.FITC.A.', y = 'X.PE.Tx.RD.A.'), subset="singlets", arrange = FALSE)
p<- p + geom_hex(bins=256)
p<- p + geom_gate(gs_get_pop_paths(auto_gs)[4:7]) 
p<- p + geom_stats(gs_get_pop_paths(auto_gs)[4:7])
p<- p + theme(strip.text = element_text(size = 7))
myPars <- ggcyto_par_set(limits = list(y = c(3,5), x = c(3,5)))
p<- p  + myPars
p

#Removing stuff
gs_pop_remove(auto_gs, "singlets")

#statistics
gs_pop_get_stats(auto_gs)
gs_pop_get_stats(auto_gs, "noneDebris_gate", "percent")
gs_pop_get_stats(auto_gs, "noneDebris_gate", type = pop.MFI)

pop.quantiles <- function(fr){
  chnls <- colnames(fr)
  res <- matrixStats::colQuantiles(exprs(fr), probs = 0.75)
  names(res) <- chnls
  res
}
gs_pop_get_stats(auto_gs, gs_get_pop_paths(auto_gs), type = pop.quantiles)

pop.mean <- function(fr){
  chnls <- colnames(fr)
  res <- colMeans(exprs(fr))
  names(res) <- chnls
  res
}
gs_pop_get_stats(auto_gs, gs_get_pop_paths(auto_gs), type = pop.mean)









#### Install the libraries #####################################################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("flowCore")
# BiocManager::install("ggplot2")
# BiocManager::install("ggpubr")
# BiocManager::install("pheatmap")
# BiocManager::install("tidyr")
# BiocManager::install("FlowRepositoryR")
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("saeyslab/FlowSOM")
devtools::install_github("saeyslab/PeacoQC")

#### Download the data #########################################################
ds <- FlowRepositoryR::download(FlowRepositoryR::flowRep.get("FR-FCM-ZZQY"), 
                                "Data/Raw")
ds <- FlowRepositoryR::download(FlowRepositoryR::flowRep.get("FR-FCM-Z2TQ"), 
                                "Data/Raw")

# If the above lines of code give an error, download the data directly from the
# FlowRepository website:
# https://flowrepository.org/experiments/833/download_ziped_files
# https://flowrepository.org/experiments/3002/download_ziped_files

#### Prepare data ##############################################################
#microbenchmark::microbenchmark({
# 1. Load the libraries
library(flowCore)
library(FlowSOM)
library(ggplot2)

# 2. Define the general and preprocessing variables
file_pattern <- "\\d.fcs" #digit at the end and fcs extension
reference_file <- read.FCS("C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15709 15Feb2022 Simone/Samples_Tube_040 NSG  in house control 1_052.fcs", truncate_max_range = FALSE)
reference_marker <- "PE-A" # Scatter values will be scaled to have the same range

markers_of_interest <- c("SSC-A", "CD8", "CD4", "CD3", "CD33",
                         "CD11b", "CD25", "mCD45", "huCD45", "CD56", 
                         "CD19", "CD66b")

live_gate <- flowCore::polygonGate(filterId = "Live",
                                   .gate = matrix(data = c(60000, 100000, 150000, 
                                                           250000, 250000, 60000, 
                                                           60000, 1.6, 1.9, 2.5,
                                                           2.5, -0.3, -0.3, 1.6),
                                                  ncol = 2,
                                                  dimnames = list(c(), 
                                                                  c("FSC-A", 
                                                                    "APC-Cy7-A"))))

# 3. Define and create the directories
dir_prepr <- "C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/Preprocessed/" #where the preprocessed data will be stored
dir_QC <- "C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/QC/" #where the data QC results will be stored
dir_RDS <- "C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/RDS/" #where the R objects will be stored
dir_results <- "C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/Results/" #where the results will be stored
dir_raw <- "C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15709 15Feb2022 Simone/" #where the raw data is located
#path_comp <- "C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/COMP/" #where comp matrix is located
comp_list <- spillover(reference_file)
compensation_matrix <- comp_list[[1]]


for (path in c(dir_prepr, dir_QC, dir_RDS, dir_results)){
  dir.create(path)
}

# 4. Prepare some additional information for preprocessing the files 
# given the variable choices of step 2.
files <- list.files(path = dir_raw,
                    pattern = file_pattern)
channels_of_interest <- GetChannels(object = reference_file,
                                    markers = markers_of_interest, 
                                    exact = FALSE)
#compensation_matrix <- read.csv(path_comp, 
#                                check.names = FALSE, row.names = 1)

colnames(compensation_matrix) <- sub(" :: .*", "",         
                                     colnames(compensation_matrix))

# Compute transformation list
ff_m <- PeacoQC::RemoveMargins(reference_file, channels_of_interest)
names(ff_m)
exprs(ff_m)
each_col(ff_m, median)

ff_c <- flowCore::compensate(ff_m, compensation_matrix)
translist <- estimateLogicle(ff_c, colnames(compensation_matrix))
ff_t <- flowCore::transform(ff_c, translist)
q5_goal <- quantile(exprs(ff_t)[,reference_marker], 0.05)
q95_goal <- quantile(exprs(ff_t)[,reference_marker], 0.95)
q5_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.05)
q95_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.95)
SSCA_a <- (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA)
SSCA_b <- q5_goal - q5_SSCA * (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA)
translist <- c(translist, 
               transformList("SSC-A", flowCore::linearTransform(a = SSCA_a,
                                                                b = SSCA_b)))

# 5. Read the first fcs file into a flowframe
ff <- read.FCS(paste0(dir_raw, files[1]), truncate_max_range = FALSE)

# 6. Remove margin events
ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest)

# 7. Compensate
ff_c <- flowCore::compensate(ff_m, compensation_matrix)

# 8. Transform, logicle for marker channels, linear for scatter channel
ff_t <- flowCore::transform(ff_c, translist)

# 9. Remove doublets and filter live cells
ff_s <- PeacoQC::RemoveDoublets(ff_t)
#selected_live <- filter(ff_s, live_gate)
#ff_l <- ff_s[selected_live@subSet, ]

# 10. QC with PeacoQC
PQC <- PeacoQC::PeacoQC(ff = ff_s,
                        channels = channels_of_interest,
                        plot = TRUE, save_fcs = FALSE,
                        output_directory = dir_QC)

# 11. Save the preprocessed data
write.FCS(PQC$FinalFF,
          file = paste0(dir_prepr, files[1]))

# 12. Visualize the preprocessing
filter_plot <- function(ff_pre, ff_post, title, channel_x, channel_y){
  df <- data.frame(x = exprs(ff_pre)[,channel_x],
                   y = exprs(ff_pre)[,channel_y])
  i <- sample(nrow(df), 10000)
  if (!"Original_ID" %in% colnames(exprs(ff_pre))) {
    ff_pre@exprs <- cbind(ff_pre@exprs,
                          Original_ID = seq_len(nrow(ff_pre@exprs)))
  }
  p <- ggplot(df[i,], aes(x = x, y = y)) +
    geom_point(size = 0.5,
               color = ifelse(exprs(ff_pre)[i,"Original_ID"] %in%
                                exprs(ff_post)[,"Original_ID"], 'blue', 'red')) +
    xlab(GetMarkers(ff_pre, channel_x)) + 
    ylab(GetMarkers(ff_pre, channel_y)) +
    theme_minimal() + theme(legend.position = "none") +
    ggtitle(title)
  return(p)
}
to_plot <- list(list(ff_pre = ff,
                     ff_post = ff_m,
                     title = "Removed margin events",
                     channel_x = "PerCP-Cy5-5-A",
                     channel_y = "BV605-A"),
                list(ff_pre = ff_t,
                     ff_post = ff_s,
                     title = "Removed doublets",
                     channel_x = "FSC-A",
                     channel_y = "FSC-H"),
                list(ff_pre = ff_s,
                     ff_post = ff_s,
                     title = "Removed debris and dead cells",
                     channel_x = "FSC-A",
                     channel_y = "APC-Cy7-A"),
                list(ff_pre = ff_l,
                     ff_post = PQC$FinalFF,
                     title = "Removed low quality events",
                     channel_x = "Time",
                     channel_y = "PerCP-Cy5-5-A"))

plot_list <- list()
for (plot in to_plot) {
  plot_list[[length(plot_list) + 1]] <- filter_plot(ff_pre = plot$ff_pre,
                                                    ff_post = plot$ff_post,
                                                    title = plot$title,
                                                    channel_x = plot$channel_x,
                                                    channel_y = plot$channel_y)
}

png(paste0(dir_QC, sub("fcs", "png", files[1])), width = 1920)
print(ggpubr::ggarrange(plotlist = plot_list, nrow = 1))
dev.off()

# 13. Run the preprocessing pipeline for all the files
for (file in files){
  ff <- read.FCS(paste0(dir_raw, file), truncate_max_range = FALSE)
  ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest)
  ff_c <- flowCore::compensate(ff_m, compensation_matrix)
  ff_t <- flowCore::transform(ff_c, translist)
  ff_s <- PeacoQC::RemoveDoublets(ff_t)
  #selected_live <- filter(ff_s, live_gate)
  #ff_l <- ff_s[selected_live@subSet, ]
  
  PQC <- PeacoQC::PeacoQC(ff = ff_s,
                          channels = channels_of_interest,
                          plot = TRUE, save_fcs = FALSE,
                          output_directory = dir_QC)
  
  write.FCS(PQC$FinalFF,
            file = paste0(dir_prepr, file))
}



###EFE DELETED PLOT CODE




# 14. Perform quality control between all files
# 14.(A) Plot the signal per channel and per file
# 14.(A)(i) Define the variables
file_names <- sub(".*15_(.*).fcs", "\\1", files)
file_groups <- rep(c("NSG", "Control"), 
                   times = c(26, 2))

# 14.(A)(ii) Make the overview plot
PlotFileScatters(input = paste0(dir_prepr, files),
                 channels = channels_of_interest,
                 names = file_names, legend = TRUE,
                 groups = file_groups, nrow = 2,
                 plotFile = paste0(dir_QC, "file_scatters.png"))

# 14.(B) Perform principal commponent analysis (PCA)
# 14.(B)(i) Retrieve the median marker expression values per file
medians <- matrix(data = NA,
                  nrow = length(files), ncol = length(channels_of_interest),
                  dimnames = list(files, channels_of_interest))

for (file in files){
  ff <- read.FCS(paste0(dir_prepr, file))
  medians[file,] <- apply(exprs(ff)[,channels_of_interest], 2, median)
}

# 14.(B)(ii) Calculate the PCs
pc <- prcomp(medians, scale. = TRUE)

# 14.(B)(iii) Visualize the PCs
ggplot(data.frame(pc$x[,1:2], file_groups)) + 
  geom_point(aes(x= PC1, y = PC2, col = file_groups)) +
  theme_minimal()
ggsave(paste0(dir_QC, "file_PCA.png"), width = 5)

#}, times = 10)

#### Create an aggregate file ##################################################
#microbenchmark::microbenchmark({

# 15. Choose the number of cells to include in the aggregate file
n <- 700000

# 16. Make an aggregate file
set.seed(2020)
agg <- AggregateFlowFrames(paste0(dir_prepr, files),
                           cTotal = n,
                           writeOutput = TRUE,
                           outputFile = paste0(dir_prepr, "aggregate.fcs"))

#}, times = 10)

#### Train FlowSOM model #######################################################
#microbenchmark::microbenchmark({

# 17. Specify the FlowSOM variables
SOM_x <- 10
SOM_y <- 10
n_meta <- 8
seed <- 2020
scaling <- FALSE


names(agg)
exprs(agg)
markers_of_interest <- c("SSC-A", "CD8", "CD4", "CD3", "CD33",
                         "CD11b", "CD25", "mCD45", "huCD45", "CD56", 
                         "CD19", "CD66b")
markers_of_interest <- c("SSC-A", "CD4", "CD3", "CD33",
                         "CD11b", "CD25", "mCD45", "huCD45", "CD56", 
                         "CD19", "CD66b")

# 18. Compute the FlowSOM object
fsom <- FlowSOM(input = agg,
                scale = scaling,
                colsToUse = c(4,7:16),
                seed = seed,
                nClus = n_meta,
                xdim = SOM_x, ydim = SOM_y)
saveRDS(fsom, paste(dir_RDS, "fsom.rds"))

# 19. Visualize the FlowSOM object
PlotStars(fsom = fsom,
          backgroundValues = fsom$metaclustering)
ggsave(paste0(dir_results, "fsom_tree.pdf"),height = 8.5, width = 11)

#}, times = 10)






PlotDimRed(fsom, colsToUse = fsom$map$colsUsed, colorBy = "metaclusters")

