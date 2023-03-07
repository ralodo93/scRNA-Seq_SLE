library(glue)
library(tidyverse)
library(Seurat)
library(harmony)
library(parallel)

merge.seurat <- readRDS("00_rds/01_filtered_cells.rds")
merge.seurat <- NormalizeData(merge.seurat, verbose = F)

merge.seurat <- CellCycleScoring(object = merge.seurat, 
                                 g2m.features = cc.genes.updated.2019$g2m.genes, 
                                 s.features = cc.genes.updated.2019$s.genes)

merge.seurat$Cycle.Score <- merge.seurat$S.Score - merge.seurat$G2M.Score

DefaultAssay(merge.seurat) <- "RNA"

merge.seurat <- FindVariableFeatures(merge.seurat)
merge.seurat <- ScaleData(merge.seurat,vars.to.regress = c("Cycle.Score"))

## Run Harmony on data
merge.seurat <- RunPCA(merge.seurat, verbose = FALSE)
merge.seurat <- RunHarmony(merge.seurat, group.by.vars = "orig.ident", reduction.save = "Harmony")

invisible(gc())
# Find significant PCs
stdv <- merge.seurat[["Harmony"]]@stdev
sum.stdv <- sum(merge.seurat[["Harmony"]]@stdev)
percent.stdv <- (stdv / sum.stdv) * 100
cumulative <- cumsum(percent.stdv)
co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                     percent.stdv[2:length(percent.stdv)]) > 0.1), 
            decreasing = T)[1] + 1
min.pc <- min(co1, co2)
print(min.pc)

merge.seurat <- RunUMAP(merge.seurat, reduction = "Harmony", dims = 1:min.pc)
merge.seurat <- FindNeighbors(merge.seurat, reduction = "Harmony", dims = 1:min.pc)
invisible(gc())

saveRDS(merge.seurat, file = "00_rds/02_dim_reduction_cells.rds")