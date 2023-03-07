library(Seurat)
library(tidyverse)
library(harmony)
library(parallel)
library(glue)
library(optparse)
library(cowplot)
library(ggrepel)

dir.create("results", showWarnings = F)

alldata  <- readRDS("00_rds/03_clustered_cells.rds")

lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

es.max = sctype_score(scRNAseqData = alldata@assays$RNA@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

invisible(gc())

cL_results = data.table::rbindlist(parallel::mclapply(unique(alldata@meta.data$integrated_snn_res.0.5), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(alldata@meta.data[alldata@meta.data$integrated_snn_res.0.5==cl, ])]), decreasing = !0)
  data.frame(cluster = cl, 
             type = names(es.max.cl), 
             scores = es.max.cl, 
             ncells = sum(alldata@meta.data$integrated_snn_res.0.5==cl))
}, mc.cores = 10)) %>% group_by(cluster) %>% mutate(zscore = (scores - mean(scores))/sd(scores))

cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores) 

assign_cluster <- c("0" = "N_CD4_T",
                    "1" = "Mono_CD14",
                    "2" = "N_CD8_T",
                    "3" = "B",
                    "4" = "CD8_NK",
                    "5" = "E_CD4_T",
                    "6" = "CD8_NK",
                    "7" = "NK",
                    "8" = "N_B",
                    "9" = "Mono_CD16",
                    "10" = "ISG_Imm",
                    "11" = "NK",
                    "12" = "Platelets",
                    "13" = "cDC",
                    "14" = "CD8_NK",
                    "15" = "Plasma",
                    "16" = "Eryth",
                    "17" = "pDC")

assign_cluster <- data.frame(cluster_number = names(assign_cluster), cluster_name = paste0("C",names(assign_cluster),"-",assign_cluster))
cL_results <- cL_results %>% left_join(assign_cluster, by = c("cluster" = "cluster_number"))
saveRDS(cL_results, file = "results/cl_scores.rds")


meta_data_add <- FetchData(alldata, c("integrated_snn_res.0.5","UMAP_1","UMAP_2")) %>% rownames_to_column("cell_id") %>% left_join(assign_cluster, by = c("integrated_snn_res.0.5" = "cluster_number")) %>% 
  column_to_rownames("cell_id")

alldata <- alldata[,meta_data_add$cell_id]
alldata@meta.data <- cbind(alldata@meta.data, meta_data_add)

# genes <- c("HBB", "PPBP","S100A8","CD14","LYZ","CST3",'CD1C',"FCER1G","FCGR3A","IL7R","CD3E","CD3D","CD8A","CD8B",
#            "GZMA","GZMB","GZMH","GZMK", "NKG7", "MS4A1", "CD79A", "MZB1","IGJ",'LILRA4',"IRF7",'IFI44L','ISG15','IFI6')
# genes <- genes[genes %in% rownames(alldata)]
# 
# expression_vals <- alldata@assays$RNA@data[genes,] %>%
#   as.data.frame() %>%
#   rownames_to_column("gene") %>%
#   pivot_longer(cols = colnames(alldata), values_to = "expression", names_to = "cell_id") %>%
#   left_join(FetchData(alldata,c("orig.ident","UMAP_1","UMAP_2",paste0("integrated_snn_res.",opt_list$resolution),"cluster_name")) %>%
#               rownames_to_column("cell_id")) %>%
#   as.data.frame()
# 
# expression_vals <- expression_vals %>%
#   mutate(cluster = as.numeric(!!as.symbol(paste0("integrated_snn_res.",opt_list$resolution)))) %>%
#   arrange(cluster) %>%
#   mutate(cluster_name = factor(cluster_name, levels = unique(cluster_name)))
# 
# plot_object <- mclapply(genes, function(id){
#   expression_vals_gene <- expression_vals %>%
#     filter(gene == id) %>%
#     as.data.frame()
#   
#   ggplot(expression_vals_gene, aes(x = UMAP_1, y = UMAP_2))+
#     geom_point(shape = ".", aes(color = expression))+
#     theme(panel.background = element_blank(),
#           axis.title = element_text(size = 5),
#           axis.text = element_blank(),
#           axis.ticks = element_blank(),
#           axis.line = element_line(size = 0.1),
#           legend.key = element_rect(fill = "white"),
#           legend.key.size = unit(0.2,"cm"),
#           legend.title = element_text(size = 5),
#           legend.text = element_text(size = 4))+
#     scale_color_gradient(low = "gray90", high = "darkred")+
#     xlab("UMAP 1")+
#     ylab("UMAP 2")+
#     labs(color = id)
# }, mc.cores = 10)
# 
# umap_plots <- plot_grid(plotlist = plot_object, ncol = 4)
# ggsave(umap_plots, filename = opt_list$umap,
#        units = "cm",height = 20, width = 20,dpi = 500,
#        bg = "white")
# 
# counts_per_cluster <- expression_vals %>%
#   select(cell_id,orig.ident,cluster_name) %>%
#   unique() %>%
#   group_by(cluster_name) %>%
#   summarise(n_cells_cluster = n()) %>%
#   as.data.frame()
# 
# dot_plot_figure <- expression_vals %>%
#   filter(expression > 0) %>%
#   group_by(cluster_name, gene) %>%
#   summarise(n_cells_pos_expression = n(),
#             mean_expression = mean(expression)) %>%
#   left_join(counts_per_cluster) %>%
#   mutate(percentage = (n_cells_pos_expression / n_cells_cluster)*100,
#          gene = factor(gene, levels = genes)) %>%
#   as.data.frame()
# 
# 
# dot_plot_matrix <- dot_plot_figure %>%
#   select(cluster_name,gene,mean_expression) %>%
#   pivot_wider(names_from = "gene", values_from = "mean_expression") %>%
#   column_to_rownames("cluster_name") %>%
#   as.data.frame()
# 
# hclust_matrix <- as.dendrogram(hclust(dist(dot_plot_matrix)))
# order_clust <- unique(dot_plot_figure$cluster_name)[order.dendrogram(hclust_matrix)]
# dot_plot_figure$clust_fact <- factor(dot_plot_figure$cluster_name, levels = order_clust)
# 
# plot_marker <- ggplot(dot_plot_figure, aes(x = gene, y = clust_fact))+
#   geom_point(aes(fill = mean_expression, size = percentage), shape = 21, color = "black", stroke = 0.2)+
#   theme(panel.background = element_rect(fill = "white"),
#         axis.title = element_blank(),
#         axis.line = element_line(),
#         legend.key = element_rect(fill = "white"),
#         axis.text.x = element_text(angle = 90, vjust = 0.5),
#         axis.ticks = element_blank())+
#   labs(fill = "Expression", size = "Percentage")+
#   scale_fill_gradient(low = "gray90", high = "darkred")
# 
# ggsave(filename = opt_list$dotplot,
#        plot = plot_marker, units = "cm",height = 20, width = 20,dpi = 500,
#        bg = "white")


message("Saving Data\n")

saveRDS(alldata, file = "00_rds/04_identification.rds")