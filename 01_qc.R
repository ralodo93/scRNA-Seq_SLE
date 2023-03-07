library(GEOquery)
library(glue)
library(tidyverse)
library(Seurat)
library(scuttle)
library(scran)
library(scater)
library(scDblFinder)
library(parallel)
library(SingleCellExperiment)
library(cowplot)

unlink("data_mtx/",recursive = T)
unlink("data/",recursive = T)
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

dir.create("data_mtx")
dir.create("data")
# download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135779&format=file","GSE135779_RAW.tar")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135779/suppl/GSE135779_genes.tsv.gz","data_mtx/features.tsv.gz")

untar("GSE135779_RAW.tar",exdir = "data_mtx")

patterns = list.files("data_mtx",pattern = "mtx.gz")
gsms = sapply(strsplit(patterns,"_"), `[`, 1)

st = getGEO("GSE135779")
phen = data.frame(phenoData(st$GSE135779_series_matrix.txt.gz)@data)

### clin_info was downloaded from https://static-content.springer.com/esm/art%3A10.1038%2Fs41590-020-0743-0/MediaObjects/41590_2020_743_MOESM3_ESM.xlsx
clin_info = read.csv("clinical_info.csv")
clin_info <- clin_info[clin_info$Groups %in% c("cSLE","cHD"),]

clin_info$SLEDAI_Group <- ifelse(clin_info$Groups == "cHD",
                                 "None", ifelse(clin_info$SLEDAI == "ND",
                                                "No Data",ifelse(as.numeric(clin_info$SLEDAI) < 4,
                                                                 "Low SLEDAI", "High SLEDAI")))

clin_info$LYMPHOCYTE_per_num <- ifelse(clin_info$Groups == "cHD",
                                       "None", ifelse(clin_info$LYMPHOCYTE_per == "ND",
                                                      "No Data",gsub(",",".",clin_info$LYMPHOCYTE_per)))
clin_info$NEUTROPHIL_per_num <- ifelse(clin_info$Groups == "cHD",
                                       "None", ifelse(clin_info$NEUTROPHIL_per == "ND",
                                                      "No Data",gsub(",",".",clin_info$NEUTROPHIL_per)))
clin_info$NLR <- ifelse(clin_info$Groups == "cHD",
                        "None", ifelse(clin_info$NEUTROPHIL_per == "No Data",
                                       "No Data",as.numeric(clin_info$NEUTROPHIL_per_num) / as.numeric(clin_info$LYMPHOCYTE_per_num)))


healthy_NLR <- read.csv("Raw_NLR_HC.csv",header = T)
n <- nrow(healthy_NLR)
healthy_NLR <- data.frame(dataset=rep("",n),
                          Cluster=rep("Healthy",n),
                          Disease=rep("",n),
                          Neu=rep("",n),
                          Lym=rep("",n),
                          NLR=as.numeric(gsub(",",".",healthy_NLR$Ratio.neutro...lympho)),
                          Interaction="Healthy")

q1 <- quantile(healthy_NLR$NLR,0.1,na.rm = T)
q3 <- quantile(healthy_NLR$NLR,0.9,na.rm = T)

clin_info$NLR_Group <- ifelse(clin_info$Groups == "cHD",
                              "None", ifelse(clin_info$NEUTROPHIL_per == "No Data",
                                             "No Data",ifelse(as.numeric(clin_info$NLR) < q1,
                                                              "Low NLR", ifelse(as.numeric(clin_info$NLR) > q3,
                                                                                "High NLR", "Normal NLR"))))

clin_info$Condition = ifelse(clin_info$Groups == "cHD","Healthy","SLE")

cc = strsplit(phen$title," \\[JB")
phen$title <- unlist(cc)[2*(1:length(phen$title))-1]
phen = merge(phen,clin_info,by.x = "title",by.y = "Names")

phen <- phen[,c("title","geo_accession","Condition","SLEDAI_Group","NLR_Group","Batch")]
gsms = phen$geo_accession
for (gsm in gsms){
  dir.create(glue("data/{gsm}"))
  files = list.files("data_mtx",gsm)
  file.copy(file.path("data_mtx", files), glue("data/{gsm}"), overwrite = TRUE)
  file.copy("data_mtx/features.tsv.gz",glue("data/{gsm}"))
  files = list.files(glue("data/{gsm}"))
  barfile = files[grep("barcodes",files)]
  file.rename(glue("data/{gsm}/{barfile}"),glue("data/{gsm}/barcodes.tsv.gz"))
  mtxfile = files[grep("mtx",files)]
  file.rename(glue("data/{gsm}/{mtxfile}"),glue("data/{gsm}/matrix.mtx.gz"))
}

write.table(phen,"data/clin_data.tsv",row.names =F, quote = F,sep = "\t")

rownames(phen) <- phen$geo_accession
datasets = phen$geo_accession
seurat_objects <- mclapply(datasets,function(dataset){
  path <- glue("data/{dataset}")
  seu.data <- Read10X(data.dir = path)
  seu <- CreateSeuratObject(counts = seu.data,
                            project = dataset,
                            min.cells = 3)
  metadata <- phen[phen$geo_accession == dataset,]
  for(metadata_column in colnames(metadata)){
    metadata_value <- metadata[,c(metadata_column)]
    seu <- AddMetaData(seu,metadata_value,col.name = metadata_column)
  }
  seu$percent.mito <- PercentageFeatureSet(seu, pattern='^MT-')
  seu$percent.ribo <- PercentageFeatureSet(seu, pattern="^RP[SL]")
  return(seu)
},mc.cores = 10)

names(seurat_objects) <- datasets
first_object <- seurat_objects[[1]]
others <- c()

for (i in 2:length(seurat_objects)){
  others <- c(others,seurat_objects[[i]])
}

alldata <- merge(first_object,others,add.cell.ids = datasets)

# metadata_00 <- alldata@meta.data %>%
#   rownames_to_column("cell_id") %>%
#   pivot_longer(cols = c(nFeature_RNA,percent.mito,percent.ribo),values_to = "value",names_to = "stats") %>% unique() %>%
#   mutate(Mantener = ifelse(stats == "nFeature_RNA" & value > 300, "Si",
#                            ifelse(stats == "percent.mito" & value < 20,"Si",
#                                   ifelse(stats == "percent.ribo" & value > 5, "Si","No")))) %>%
#   mutate(stats = ifelse(stats == "nFeature_RNA", "Número de genes expresados",
#                         ifelse(stats == "percent.mito","Porcentaje de genes mitocondriales", "Porcentaje de genes ribosomales")))
# 
# c1 <- ggplot(metadata_00 %>% filter(stats == "Número de genes expresados"), aes(x = orig.ident, y = value, color = Mantener))+
#   theme_void()+
#   theme(legend.position = "none",
#         axis.title.y = element_text(size = 6, angle = 90, hjust = 1, margin = margin(r = 2)),
#         axis.text.y = element_text(size = 5),
#         plot.margin = margin(5,5,5,5),
#         axis.line = element_line(linewidth = 0.2))+
#   geom_jitter(size = 0.01)+
#   ylab("Número de genes expresados")+
#   scale_color_manual(values = c("Si" = "#4caf50", "No" = "#f44336"))+
#   scale_y_continuous(expand = c(0,0))
# 
# c2 <- ggplot(metadata_00 %>% filter(stats == "Porcentaje de genes mitocondriales"), aes(x = orig.ident, y = value, color = Mantener))+
#   theme_void()+
#   theme(legend.position = "none",
#         axis.title.y = element_text(size = 6, angle = 90, hjust = 1, margin = margin(r = 3)),
#         axis.text.y = element_text(size = 5),
#         plot.margin = margin(5,5,5,5),
#         axis.line = element_line(linewidth = 0.2))+
#   geom_jitter(size = 0.01)+
#   ylab("Porcentaje de genes mitocondriales")+
#   scale_color_manual(values = c("Si" = "#4caf50", "No" = "#f44336"))+
#   scale_y_continuous(expand = c(0,0))
# 
# c3 <- ggplot(metadata_00 %>% filter(stats == "Porcentaje de genes ribosomales"), aes(x = orig.ident, y = value, color = Mantener))+
#   theme_void()+
#   theme(axis.title.y = element_text(size = 6, angle = 90, hjust = 1, margin = margin(r = 3)),
#         axis.text.y = element_text(size = 5),
#         axis.title.x = element_text(size = 6, hjust = 1, margin = margin(t = 3)),
#         plot.margin = margin(5,5,5,5),
#         axis.line = element_line(linewidth = 0.2))+
#   geom_jitter(size = 0.01)+
#   ylab("Porcentaje de genes ribosomales")+
#   xlab("Individuos")+
#   scale_color_manual(values = c("Si" = "#4caf50", "No" = "#f44336"))+
#   scale_y_continuous(expand = c(0,0))
# 
# legend <- get_legend(
#   # create some space to the left of the legend
#   c3 + theme(legend.box.margin = margin(0, 0, 0, 12))
# )
# 
# prow <- egg::ggarrange(c1, c2, c3 + theme(legend.position = "none"),
#                        nrow = 1)

# meta_data <- FetchData(alldata,c("orig.ident","nFeature_RNA","nCount_RNA","percent.mito","percent.ribo"))

# raw_data <- mclapply(c("nFeature_RNA","nCount_RNA","percent.mito","percent.ribo"), function(meta_col){
#   meta_data_col <- meta_data %>%
#     select(orig.ident, all_of(meta_col))
#   colnames(meta_data_col)[2] <- "meta_col"
#   ggplot(meta_data_col, aes(x = orig.ident, y = meta_col))+
#     geom_jitter(size = 0.2, shape = 21, color = "black", fill = "deepskyblue", stroke = 0.05)+
#     theme(axis.text.x = element_blank(),
#           axis.text.y = element_text(size = 5),
#           axis.title = element_text(size = 6),
#           panel.background = element_blank(),
#           axis.line = element_line())+
#     xlab("Samples")+
#     ylab(meta_col)
# })
# 
# raw_data <- plot_grid(plotlist = raw_data)
# 
# ggsave(filename =  "results/01_qc_before_filtering.png",
#        plot = raw_data, units = "cm",height = 15, width = 15,dpi = 500,
#        bg = "white")

alldata<- subset(alldata, nFeature_RNA > 300)
alldata<- subset(alldata, percent.mito < 20)
alldata<- subset(alldata, percent.ribo > 0.05)

alldata<- alldata[!grepl("MALAT1",rownames(alldata)),]

dta <- as.SingleCellExperiment(alldata)
dta.cnt <- logNormCounts(dta)
dec <- modelGeneVar(dta.cnt, block = dta.cnt$orig.ident) ## este seria el id de muestra
hvgs = getTopHVGs(dec, n = 2000)
dta.cnt <- runPCA(dta.cnt, subset_row = hvgs)
dta.cnt <- runUMAP(dta.cnt, pca = 10)
dbl.dens <- computeDoubletDensity(dta.cnt,d=ncol(reducedDim(dta.cnt)))
summary(dbl.dens)
#dta.cnt$DoubletScore <- dbl.dens
#plotUMAP(dta.cnt, colour_by="DoubletScore",point_size=0.01)
dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),
                                 method="griffiths", returnType="call")
summary(dbl.calls)
alldata@meta.data$doublets<-dbl.calls

# meta_data <- FetchData(alldata, c("orig.ident","doublets"))

# geom_count <- ggplot(meta_data, aes(x = orig.ident, fill = doublets))+
#   geom_bar(color = "black", position = "fill")+
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 8),
#         axis.title = element_text(size = 9),
#         panel.background = element_blank(),
#         axis.line = element_line())+
#   scale_y_continuous(expand = c(0,0))+
#   xlab("Samples")+
#   ylab("Proportion of singlet and doublets\nper sample")+
#   scale_fill_brewer(type = "qual")+
#   labs(fill = "")
# 
# ggsave(filename =  "results/02_singlet_and_doublets.png",
#        plot = geom_count, units = "cm",height = 15, width = 15,dpi = 500,
#        bg = "white")

alldata<- subset(alldata, doublets == "singlet")

# meta_data <- FetchData(alldata, c("orig.ident"))
# geom_count <- ggplot(meta_data, aes(x = orig.ident))+
#   geom_bar(fill = "steelblue", color = "black")+
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 8),
#         axis.title = element_text(size = 9),
#         panel.background = element_blank(),
#         axis.line = element_line())+
#   scale_y_continuous(expand = c(0,0))+
#   xlab("Samples")+
#   ylab("Numer of cells\nafter quality control")
# 
# ggsave(filename =  "results/03_n_cells_after_qc.png",
#        plot = geom_count, units = "cm",height = 15, width = 15,dpi = 500,
#        bg = "white")

# metadata_01 <- alldata@meta.data %>%
#   select(orig.ident, Condition) %>%
#   mutate(Condition = ifelse(Condition == "SLE", "LES","Sanos")) %>%
#   group_by(orig.ident, Condition) %>%
#   summarise(n = n())
# 
# c2 <- ggplot(metadata_01, aes(x = Condition, y = n, fill = Condition))+
#   theme_void()+
#   theme(axis.title.y = element_text(size = 6, angle = 90, hjust = 1, margin = margin(r = 3)),
#         axis.text.y = element_text(size = 5),
#         axis.title.x = element_text(size = 6, hjust = 1, margin = margin(t = 3)),
#         axis.text.x = element_text(size = 6, margin = margin(t = 3)),
#         plot.margin = margin(5,5,5,5),
#         axis.line = element_line(linewidth = 0.2),
#         legend.position = "none")+
#   geom_jitter(size = 2, shape = 21, color = "black", stroke = 0.1)+
#   ylab("Número de células por muestra")+
#   xlab("Condición")+
#   scale_fill_manual(values = c("Sanos" = "#7eb0d5", "LES" = "#b2e061"))
# 
# metadata_01 <- metadata_01 %>% group_by(Condition) %>% summarise(sum = sum(n))
# 
# c3 <- ggplot(metadata_01, aes(x = Condition, y = sum, fill = Condition))+
#   geom_col(color = "black", linewidth = 0.2)+
#   ylab("Número de células totales")+
#   xlab("Condición")+
#   theme_void()+
#   theme(axis.title.y = element_text(size = 6, angle = 90, hjust = 1, margin = margin(r = 3)),
#         axis.text.y = element_text(size = 5),
#         axis.title.x = element_text(size = 6, hjust = 1, margin = margin(t = 3)),
#         axis.text.x = element_text(size = 6, margin = margin(t = 3)),
#         plot.margin = margin(5,5,5,5),
#         axis.line = element_line(linewidth = 0.2),
#         legend.position = "none")+
#   scale_fill_manual(values = c("Sanos" = "#7eb0d5", "LES" = "#b2e061"))+
#   scale_y_continuous(expand = c(0,0))
# 
# fig1 <- ggdraw()+
#   draw_plot(prow,x = 0.01, y = 0.5, height = 0.5, width = 0.99)+
#   draw_plot(c2, x = 0.01, y = 0, height = 0.5, width = 0.49)+
#   draw_plot(c3, x = 0.51, y = 0, height = 0.5, width = 0.49)+
#   draw_plot_label(c("a","b","c"),x = c(0,0,0.5), y = c(1,0.5,0.5), size = 8,fontface = "plain")
# 
# ggsave(plot = fig1,filename = "Figura Tesis 1.tiff", height = 10, width = 16.5, units = "cm", dpi=800, compression = "lzw+p", bg = "white")


dir.create("00_rds",showWarnings = F)
message("Saving Seurat Data")
print(dim(alldata))
saveRDS(alldata, file = "00_rds/01_filtered_cells.rds")