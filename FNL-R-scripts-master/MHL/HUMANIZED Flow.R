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
library(flowCore)
library(FlowSOM)
library(ggplot2)

########### 1. Generate Counts From Manual Gating ########### 
########### 1. Generate Counts From Manual Gating ###########
########### 1. Generate Counts From Manual Gating ###########
data_dir <- "C:/Users/edmondsonef/Desktop/Humanized/Flow/"
study_dir1 <- "1-05Jan2022/"
study_dir2 <- "2-02Feb2022/"
study_dir3 <- "3-15Feb2022/"
study_dir4 <- "4-28Feb2022/"
study_dir5 <- "5-02Mar2022/"
study_dir6 <- "6-10Mar2022/"
ws <- open_flowjo_xml(paste0(data_dir,study_dir,"15719 02Mar2022 Simone.wsp"))
ws
fj_ws_get_samples(ws, group_id = 4)
gs <- flowjo_to_gatingset(ws, name = 4, path=dir)
gh_pop_get_stats(gs[[1]], "/scatter/sing")
cs <- gs_pop_get_data(gs, "/scatter/sing")
fs <- cytoset_to_flowSet(cs)
gs <- flowjo_to_gatingset(wsp, name = 1, path ="C:/Users/edmondsonef/Desktop/Humanized/Flow/3-15Feb2022/")
plot(gs)
autoplot(gs[[1]])

gs_get_pop_paths(gs)
recompute(gs)

gs_pop_get_gate(gs, "/scatter")
x <- gs_pop_get_gate(gs[[1]], "/scatter")


sampStats <- gs_pop_get_stats(gs)

### Get Stats on Manual Gates or Nodes
### Get Stats on Manual Gates or Nodes
### Get Stats on Manual Gates or Nodes
### Get Stats on Manual Gates or Nodes
nodes <- c("/scatter/sing/", "/scatter/sing/hCD45+",
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
           "Q35: CD33+ , CD11b [AF647]-",
           "Q39: CD25+ , CD3-")
gs_pop_get_stats(gs, nodes, "percent")
nodeCount <- gs_pop_get_stats(gs, nodes, "count")
nodeCount
write.csv(nodeCount, "C:/Users/edmondsonef/Desktop/nodeCount.csv")



# 2. Define the general and preprocessing variables
data_dir <- "C:/Users/edmondsonef/Desktop/FACS/FLOWSET/"
file_pattern <- "\\d.fcs" #digit at the end and fcs extension
reference_file <- read.FCS(paste0(data_dir, 'Samples_Tube_015 Animal 142_027.fcs'), truncate_max_range = FALSE)
reference_marker <- "PE-A" # Scatter values will be scaled to have the same range

markers_of_interest <- c("BB515-A",
                         "BB700-P-A",
                         "APC-A",
                         "APC-Cy7-A",
                         "BV421-A",
                         "BV786-A",
                         "BUV395-A",
                         "BUV805-A",
                         "PE-A",
                         "PE-CF594-A",
                         "PE-Cy7-A")

# 3. Define and create the directories
dir_prepr <- "C:/Users/edmondsonef/Desktop/FACS/FLOWSET/1-Preprocessed/" #where the preprocessed data will be stored
dir_QC <- "C:/Users/edmondsonef/Desktop/FACS/FLOWSET/2-QC/" #where the data QC results will be stored
dir_RDS <- "C:/Users/edmondsonef/Desktop/FACS/FLOWSET/3-RDS/" #where the R objects will be stored
dir_results <- "C:/Users/edmondsonef/Desktop/FACS/FLOWSET/4-Results/" #where the results will be stored
dir_raw <- "C:/Users/edmondsonef/Desktop/FACS/FLOWSET/" #where the raw data is located
path_comp <- "C:/Users/edmondsonef/Desktop/0-Autospill/table_compensation/autospill_compensation.csv" #where comp matrix is located

for (path in c(dir_prepr, dir_QC, dir_RDS, dir_results)){
  dir.create(path)
}

# 4. Prepare some additional information for preprocessing the files 
# given the variable choices of step 2.
files <- list.files(path = dir_raw, pattern = "Samples")
files
channels_of_interest <- GetChannels(object = reference_file,
                                    markers = markers_of_interest, 
                                    exact = FALSE)
compensation_matrix <- read.csv(path_comp, 
                                check.names = FALSE, row.names = 1)

colnames(compensation_matrix) <- sub(" :: .*", "",         
                                     colnames(compensation_matrix))

# Compute transformation list
ff_m <- PeacoQC::RemoveMargins(reference_file, channels_of_interest)
names(ff_m)
#exprs(ff_m)
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
translist
#rm(ff_c, ff_m, ff_t, reference_file)

# 5. Read the first fcs file into a flowframe
ff <- read.FCS(paste0(dir_raw, files[2]), truncate_max_range = FALSE)

# 6. Remove margin events
ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest)

# 7. Compensate
ff_c <- flowCore::compensate(ff_m, compensation_matrix)

# 8. Transform, logicle for marker channels, linear for scatter channel
ff_t <- flowCore::transform(ff_c, translist)

# 9. Remove doublets and filter live cells
ff_s <- PeacoQC::RemoveDoublets(ff_t)
#selected_live <- flowCore::filter(ff_s, live_gate)
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
                     channel_y = "FSC-H"))

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

