library(Seurat)
library(tidyverse)
library(harmony)
library(parallel)
library(glue)
library(optparse)

dir.create("results", showWarnings = F)

option_list <- list(
  make_option("--input",type = "character", 
              default = "00_rds/02_dim_reduction_cells.rds", help = "Select file to reduction"),
  make_option("--output",type = "character",
              default = "00_rds/03_clustered_cells.rds", help = "Select file to out reduction"),
  make_option("--umap",type = "character",
              default = "results/04_clusters.png", help = "Select file to out reduction"))

opt_list <- parse_args(OptionParser(option_list=option_list))

data_seurat <- readRDS("00_rds/02_dim_reduction_cells.rds")

resolutions <- c(0.01,0.05,seq(0.1,1,0.1))

cat("\nFind clusters for some resolutions\n")
library(data.table)
data_cluster <- rbindlist(mclapply(resolutions, function(resolution){
  print(resolution)
  cluster <- FindClusters(object = data_seurat,resolution = resolution)
  data_cluster <- FetchData(cluster,c("seurat_clusters")) %>%
    rownames_to_column("cell_id") %>%
    mutate(resolution = resolution)
  rm(cluster)
  invisible(gc())
  return(data_cluster)
}, mc.cores = 10))

for (resolution_id in resolutions){
  insert_meta <- data_cluster[data_cluster$resolution == resolution_id,]
  rownames(insert_meta) <- insert_meta$cell_id
  insert_meta <- insert_meta[,"seurat_clusters",drop = F]
  colnames(insert_meta) <- glue("integrated_snn_res.{resolution_id}")
  data_seurat@meta.data <- cbind(data_seurat@meta.data,insert_meta)
}

# library(cowplot)
# 
# resolution_cols <- colnames(data_seurat@meta.data)[grepl("integrated_snn_res", colnames(data_seurat@meta.data))]
# resolution_vals <- gsub("integrated_snn_res.","", resolution_cols)
# 
# plot_list <- lapply(1:length(resolution_cols), function(i){
#   resolution_col <- resolution_cols[i]
#   resolution_val <- resolution_vals[i]
#   
#   meta_data <- FetchData(data_seurat, c(resolution_col, "UMAP_1","UMAP_2")) %>%
#     rename("seurat_clusters" = resolution_col)
#   
#   centroids <- meta_data %>%
#     group_by(seurat_clusters) %>%
#     summarise(meanX = median(UMAP_1),
#               meanY = median(UMAP_2)) %>%
#     as.data.frame()
#   
#   ggplot(meta_data, aes(x = UMAP_1, y = UMAP_2))+
#     geom_point(shape = ".", aes(color = seurat_clusters))+
#     theme(panel.background = element_blank(),
#           axis.title = element_blank(),
#           axis.text = element_blank(),
#           axis.ticks = element_blank(),
#           axis.line = element_line(),
#           legend.position = "none",
#           plot.title = element_text(hjust = 0.5))+
#     geom_text(data = centroids, aes(x = meanX, y = meanY, label = seurat_clusters), size = 2)+
#     ggtitle(resolution_val)
# })
# 
# # Viasualization of clusters
# plot_list <- plot_grid(plotlist = plot_list)
# 
# ggsave(plot_list, filename = opt_list$umap,
#        units = "cm",height = 20, width = 20,dpi = 500,
#        bg = "white")

cat("\nSave data_seurat\n")
saveRDS(data_seurat, file = "00_rds/03_clustered_cells.rds")