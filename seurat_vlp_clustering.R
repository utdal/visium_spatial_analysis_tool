library(Seurat)
library(ggplot2)
library(dplyr)
library(future)
plan(sequential)

#######################
# violin plots function
#######################
vln_plot_function <- function(Vln_plot_data_dir_path, project_name) {
  data_list <- list()
  file_list <- sort(list.files(path=Vln_plot_data_dir_path))
  denotion_list <- list("All Neurons", "Other", "Single", "Surrounding")
  
  for (i in seq_along(file_list)) {
    file_path <- file.path(Vln_plot_data_dir_path, file_list[i])
    data <- read.csv(file_path)
    data$Gene.Name <- make.unique(data$Gene.Name, sep = "_dup")
    rownames(data) <- data$Gene.Name
    data$Gene.Name <- NULL
    seurat_data_obj <- CreateSeuratObject(counts = data, 
                                          project = denotion_list[[i]])
                         # denotion_list[[which(i %in% denotion_list)]]
    data_list[[i]] <- seurat_data_obj
  }
  
  object.combined <- merge(data_list[[1]], 
                           y = c(data_list[[2]], data_list[[3]], data_list[[4]]), 
                           add.cell.ids = c("Multiple", "Other", "Single", "Surrounding"), 
                           project = project_name)
  
  VlnPlot(object.combined, 
          features = c("nFeature_RNA", "nCount_RNA"), 
          group.by = 'orig.ident') 
}

vln_plot_function("C:/Users/NXI220005/Desktop/Processed_files/o2", "Pig_DRG")


############
# clustering
############

clustering_function <- function(metadata_filepath) {
  clust_list = list()
  metadata_df <- read.csv(metadata_filepath)
  
  dir_path <- file.path(dirname(metadata_filepath))
  file_name <- basename(metadata_filepath)
  plot_path <- paste(dir_path, "/Seurat_processed_files", sep = "")
  # Creating a new directory if "Seurat_processed_files" dir is not already 
  # created at the given folder level where the metadata file is present
  
  if (!file.exists(plot_path)) {
    dir.create(plot_path)
  }
  
  for (i in 1:nrow(metadata_df)) {
    single_neur_fp <- metadata_df[i, "Input_to_Cluster"]
    data <- read.csv(single_neur_fp, row.names=1)
    data <- as.data.frame(data)
    s <- CreateSeuratObject(counts = data, project = metadata_df[i, "Sample_ID"],
                            min.cells = 3, min.features = 200)
    s$sex <- metadata_df[i, "Sex"]
    s$round <- metadata_df[i, "Batch"]
    s$id <- metadata_df[i, "Sample_ID"]
    s$replicate <- metadata_df[i, "Replicate"]
    s$condition <- metadata_df[i, "Condition"]
    s <- NormalizeData(s, verbose = FALSE)
    s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 2000)
    clust_list <- c(clust_list, s)
  }
  clust_list
  # Perform integration
  drg.anchors <- FindIntegrationAnchors(object.list = clust_list, dims = 1:10)
  drg.combined <- IntegrateData(anchorset = drg.anchors, dims = 1:10)
  
  # QC by sample-id
  vln_plotpath <- paste(plot_path, "/vln_plot_drg_combined.pdf", sep = "")
  vlnplot <- VlnPlot(drg.combined, features = c("nFeature_RNA", "nCount_RNA"), pt.size=0.1)
  ggsave(vln_plotpath, vlnplot)
  vln_plotpath
  
  scatter_plotpath <- paste(plot_path, "/scatter_plot_drg_combined.pdf", sep = "")
  featurescatter <- FeatureScatter(drg.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave(scatter_plotpath, featurescatter)
  
  # Perform an integrated analysis
  DefaultAssay(drg.combined) <- "integrated"
  
  # Run the standard workflow for visualization and clustering
  drg.combined <- ScaleData(drg.combined, verbose = FALSE)
  drg.combined <- RunPCA(drg.combined, 
                         features = VariableFeatures(object = drg.combined), 
                         npcs = 50)
  
  dim_plotpath <- paste(plot_path, "/dim_plot_drg_combined.pdf", sep = "")
  dimplot <- DimPlot(drg.combined, reduction = "pca")
  ggsave(dim_plotpath, dimplot)
  dim_heatmap_plotpath <- paste(plot_path, "/dim_heatmap_plot_drg_combined.pdf", sep = "")
  dimheatmap <- DimHeatmap(drg.combined, dims = 1:10, cells = 500, balanced = TRUE)
  ggsave(dim_heatmap_plotpath, dimheatmap)
  
  # UMAP and Clustering
  drg.combined <- RunUMAP(drg.combined, reduction = "pca", dims = 1:10)
  drg.combined <- FindNeighbors(drg.combined, reduction = "pca", dims = 1:10)
  drg.combined <- FindClusters(drg.combined, resolution = 1)
  umap_dim_plotpath <- paste(plot_path, "/umap_dim_plot_drg_combined.pdf", sep = "")
  umapdimplot <- DimPlot(drg.combined, reduction = "umap", label = TRUE, pt.size=1)
  ggsave(umap_dim_plotpath, umapdimplot)
  
  # at this point, you should look at markers for each cluster and at the expression
  # of genes of interest. If needed re-run clustering and change parameters where needed.
  
  # save seurat object
  seurat_filepath <- paste(plot_path, "/seurat_object_snapshot.rds", sep = "")
  saveRDS(drg.combined, file = seurat_filepath)
  
  # Get markers for each cluster
  drg.combined.markers <- FindAllMarkers(drg.combined, only.pos = TRUE, assay='RNA',
                                         min.pct = 0.25, logfc.threshold = 0.25)
  drg.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
  drg_marker_savepath <- paste(plot_path, "/drg_combined_markers_per_cluster_labeled.csv", sep = "")
  write.csv(drg.combined.markers, drg_marker_savepath)
}


###############################################################################
###############################################################################
###############################################################################
# Looking at genes
DefaultAssay(drg.combined) <- "RNA"
feature_plotpath <- paste(plot_path, "/feature_plot_drg_combined.pdf", sep = "")
featureplot <- FeaturePlot(drg.combined, features = ('SCN9A'))
ggsave(feature_plotpath, featureplot)

# Get average per cluster
average_cluster = AverageExpression(object = drg.combined, assays="RNA")
drg_avg_savepath <- paste(plot_path, "/drg_avg_cluster.csv", sep = "")
write.csv(average_cluster, drg_avg_savepath)

# Get all counts 
df6 <- GetAssayData(object = drg.combined, slot = "data", assay= "RNA")
drg_norm_count_savepath <- paste(plot_path, "/seurat_norm_counts.csv", sep = "")
write.csv(df6, drg_norm_count_savepath)

prop.table(table(Idents(drg.combined)))
WhichCells(drg.combined)
d1=FetchData(drg.combined, vars='ident')
write.csv(d1, "barcode_clusterID.csv")


###############################################################################
###############################################################################
###############################################################################

clustering_function("C:/Users/NXI220005/Desktop/vis_metadata.csv")




























#########################################
#########################################
#########################################

clustering_data_dir_path <- "C:/Users/NXI220005/Desktop/visium_results/Processed_files/Seurat_input_to_cluster"
single_neur_filelist <- sort(list.files(clustering_data_dir_path))

for (i in seq_along(single_neur_filelist)) {
  data_filepath <- file.path(clustering_data_dir_path, single_neur_filelist[i])
  print(file.path(clustering_data_dir_path, single_neur_filelist[i]))
  
  data <- read.csv(data_filepath, row.names=1)
  data <- as.data.frame(data)
  s <- CreateSeuratObject(counts = data, project = "sample1", min.cells = 3, min.features = 200)
  s$sex <- "FEMALE"     # sex
  s$round <- "ROUND1"   # sequencing run
  s$id <- "s1r1"
  s <- NormalizeData(s, verbose = FALSE)
  s <- FindVariableFeatures(s, selection.method = "vst", nfeatures = 2000)
}
#########################################
#########################################
#########################################
# Read the CSV file into a data frame
data <- read.csv("path/to/file.csv")

# Iterate over each row in the data frame
for (i in 1:nrow(data)) {
  # Access the values in each column for the current row
  col1_value <- data[i, "col1"]
  col2_value <- data[i, "col2"]
  
  # Do something with the values
  print(paste("Row", i, "values:", col1_value, col2_value))
}
#########################################
#########################################
#########################################
# Get the operating system name
os_name <- Sys.info()["sysname"]

# Construct the file path using the appropriate directory separator
if (os_name == "Windows") {
  file_path <- paste("C:", "Users", "username", "Documents", "file.txt", sep = "\\")
} else {
  file_path <- file.path("home", "username", "Documents", "file.txt")
}

# Print the file path
print(file_path)

