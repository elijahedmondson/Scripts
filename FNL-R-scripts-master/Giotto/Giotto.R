library(Giotto)

instrs = createGiottoInstructions(show_plot = FALSE,
                                  save_plot = TRUE,
                                  save_dir = 'giotto_results',
                                  python_path = "/usr/local/bin/python3.8")

qp = read.delim("qupath_quants.tsv.gz")

qp_expr = qp[,grepl("\\.\\.(Cell|Nucleus|Cytoplasm|Membrane)..",names(qp)) & !grepl("Autofluorescence|DAPI",names(qp))]
qp_expr = t(qp_expr)
colnames(qp_expr) = rownames(qp)

qp_spatial_loc = qp[,c("Centroid.X.µm","Centroid.Y.µm")]
qp_spatial_loc$Centroid.Y.µm = - qp_spatial_loc$Centroid.Y.µm
qp_spatial_loc$cell_ID = rownames(qp)
qp_spatial_loc = qp_spatial_loc[,c(3,1,2)]

gobj <- createGiottoObject(expression = qp_expr,
                           spatial_locs = qp_spatial_loc,
                           instructions = instrs)

# optionally add QuPath metadata such as marker positivity
qp_metadata = qp[,grepl("Class|phenotype",names(qp))]
qp_metadata$cell_ID = rownames(qp)

gobj<-addCellMetadata(gobj, new_metadata = qp_metadata,
                      by_column = T,
                      column_cell_ID = "cell_ID")

# Proceed with Giotto analysis such as clustering and spatial correlations
...

#Once done with the Giotto clustering analysis its results can be saved to a tsv/csv file and fed back to QuPath.


for_qp = qp[,c("Centroid.X.µm","Centroid.Y.µm")]
for_qp$Image = "giotto"
for_qp$leiden = gobj@cell_metadata$cell$rna$leiden
write.table(for_qp,"qp_giotto.tsv",quote = FALSE,sep = "\t",row.names = FALSE)