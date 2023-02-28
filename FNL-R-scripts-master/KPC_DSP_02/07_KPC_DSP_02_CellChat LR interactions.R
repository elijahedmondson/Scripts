devtools::install_github("sqjin/CellChat")
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)


load("C:/Users/edmondsonef/Desktop/DSP GeoMx/Results/KPC_seurat_2023.RData")

# Prepare input data for CelChat analysis
data.input = GetAssayData(visium.brain, slot = "data", assay = "SCT") # normalized data matrix
meta = data.frame(labels = Idents(visium.brain), row.names = names(Idents(visium.brain))) # manually create a dataframe consisting of the cell labels
unique(meta$labels) # check the cell labels
#> [1] L2/3 IT Astro   L6b     L5 IT   L6 IT   L6 CT   L4      Oligo  
#> Levels: Astro L2/3 IT L4 L5 IT L6 IT L6 CT L6b Oligo


#####
##### Create a CellChat object


# USERS can create a new CellChat object from a data matrix or Seurat. If input is a Seurat object, the meta data in the object will be used by default and USER must provide group.by to define the cell groups. e.g, group.by = "ident" for the default cell identities in Seurat object.
# 
# NB: If USERS load previously calculated CellChat object (version < 1.6.0), please update the object via updateCellChat


# load spatial imaging information
# Spatial locations of spots from full (NOT high/low) resolution images are required
spatial.locs = GetTissueCoordinates(visium.brain, scale = NULL, cols = c("imagerow", "imagecol")) 
# Scale factors and spot diameters of the full resolution images 
scale.factors = jsonlite::fromJSON(txt = file.path("/Users/jinsuoqin/Mirror/CellChat/tutorial/spatial_imaging_data_visium-brain", 'scalefactors_json.json'))
scale.factors = list(spot.diameter = 65, spot = scale.factors$spot_diameter_fullres, # these two information are required
                     fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, lowres = scale.factors$tissue_lowres_scalef # these three information are not required
)
# USER can also extract scale factors from a Seurat object, but the `spot` value here is different from the one in Seurat. Thus, USER still needs to get the `spot` value from the json file. 

###### Applying to different types of spatial imaging data ######
# `spot.diameter` is dependent on spatial imaging technologies and `spot` is dependent on specific datasets


#####
#####  Set the ligand-receptor interaction database

# Our database CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60% of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions, including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions and 16.5% of cell-cell contact interactions.
# 
# Users can update CellChatDB by adding their own curated ligand-receptor pairs.Please check our tutorial on how to do it.


cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, scale.factors = scale.factors)
#> [1] "Create a CellChat object from a data matrix"
#> Create a CellChat object from spatial imaging data... 
#> Set cell identities for the new CellChat object 
#> The cell groups used for CellChat analysis are  Astro L2/3 IT L4 L5 IT L6 IT L6 CT L6b Oligo
cellchat
#> An object of class CellChat created from a single dataset 
#>  648 genes.
#>  1073 cells. 
#> CellChat analysis of spatial data! The input spatial locations are 
#>                    imagerow imagecol
#> AAACAGAGCGACTCCT-1     3164     7950
#> AAACCGGGTAGGTACC-1     6517     3407
#> AAACCGTTCGTCCAGG-1     7715     4371
#> AAACTCGTGATATAAG-1     4242     9258
#> AAAGGGATGTAGCAAG-1     4362     5747
#> AAATAACCATACGGGA-1     3164     7537


#####
#####  Set the ligand-receptor interaction database
# Our database CellChatDB is a manually curated database of literature-supported ligand-receptor interactions in both human and mouse. CellChatDB in mouse contains 2,021 validated molecular interactions, including 60% of secrete autocrine/paracrine signaling interactions, 21% of extracellular matrix (ECM)-receptor interactions and 19% of cell-cell contact interactions. CellChatDB in human contains 1,939 validated molecular interactions, including 61.8% of paracrine/autocrine signaling interactions, 21.7% of extracellular matrix (ECM)-receptor interactions and 16.5% of cell-cell contact interactions.
# 
# Users can update CellChatDB by adding their own curated ligand-receptor pairs.Please check our tutorial on how to do it.

CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use










