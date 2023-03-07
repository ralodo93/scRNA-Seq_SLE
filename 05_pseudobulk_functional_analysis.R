library(Seurat)

data_seurat <- readRDS("00_rds/04_identification.rds")

get_pseudobulk <- function(data_seurat, sample_col = "orig.ident", groups_col = "integrated_snn_res.0.5"){
  require(tidyverse)
  require(parallel)
  require(Seurat)
  require(glue)
  metadata_cols <- data_seurat@meta.data %>% select(!!c(sample_col,groups_col))
  samples <- unique(metadata_cols %>% pull(!!sample_col))
  groups <- unique(metadata_cols %>% pull(!!groups_col))
  nrow_samples_groups <- length(samples) * length(groups)
  message("Preparing Pseudobulk")
  pseudo_bulk <- do.call("cbind", mclapply(samples, function(sample){
    sample_groups <- lapply(groups, function(group){
      index <- glue("{sample}__{group}")
      cells <- metadata_cols %>% filter(!!as.symbol(sample_col) == sample, !!as.symbol(groups_col) == group) %>% rownames_to_column("cell_id") %>% pull(cell_id)
      if (length(cells) < 10){
        return(NA)
      } else{
        tmp <- data_seurat[,cells]@assays$RNA@counts
        prop <- rowSums(tmp != 0) / nrow(tmp)
        profile <- rowSums(tmp)
        res <- data.frame(index = profile)
        colnames(res)[1] <- index
        return(res) 
      }
    })
    sample_groups <- do.call("cbind",Filter(function(a) any(!is.na(a)), sample_groups))
  }, mc.cores = 10))
  
  pseudo_bulk_seurat <- CreateSeuratObject(counts = pseudo_bulk)
  
  col.name = "Condition"
  unique_cols <- unlist(lapply(colnames(data_seurat@meta.data), function(col.name){
    if (col.name == sample_col){
      return(NA)
    } else{
      unique_rows <- data_seurat@meta.data %>% select(!!as.symbol(sample_col),col.name) %>% unique()
      if (nrow(unique_rows) == length(samples)){
        return(col.name)
      } else{
        return(NA)
      }
    }
  })) %>% na.omit()
  
  meta_data <- pseudo_bulk_seurat@meta.data %>% rownames_to_column("sample_group") %>% 
    mutate(!!sym(groups_col) := sapply(strsplit(sample_group,"__"), `[`, 2)) %>%
    left_join(data_seurat@meta.data %>% select(!!as.symbol(sample_col),unique_cols) %>% unique())
  pseudo_bulk_seurat@meta.data <- meta_data %>% column_to_rownames("sample_group")
  return(pseudo_bulk_seurat)
}

perform_functional_analysis <- function(pseudo_bulk_seurat_in, sample_col = "orig.ident",
                                        groups_col = "integrated_snn_res.0.5", comparison = "Condition",
                                        ident1 = "SLE", rem = F){
  require(dorothea)
  require(DESeq2)
  require(progeny)
  require(tidyverse)
  require(parallel)
  require(Seurat)
  require(glue)
  require(data.table)
  
  
  
  
  hs_regulon <- dorothea_hs %>%
    filter(confidence <= "C", target %in% rownames(data_seurat@assays$RNA@counts)) %>%
    as.data.frame()
  hs_regulon_list <- hs_regulon %>% df2regulon()
  
  
  groups <- pseudo_bulk_seurat_in@meta.data %>% pull(!!as.symbol(groups_col)) %>% unique()
  if (!rem){
    print("Differential expression analysis between two conditions")
    ident2 <- pseudo_bulk_seurat_in@meta.data %>% pull(!!as.symbol(comparison)) %>% unique()
    ident2 <- ident2[ident2 != ident1]
    diff_expression <- mclapply(groups, function(group){
      samples <- pseudo_bulk_seurat_in@meta.data %>% filter(!!as.symbol(groups_col) == group) %>% rownames_to_column("sample_id") %>% pull(sample_id)
      pseudo_bulk_seurat_group <- pseudo_bulk_seurat_in[,samples]
      cells1 <- nrow(pseudo_bulk_seurat_group@meta.data %>% filter(!!as.symbol(comparison) == ident1))
      cells2 <- nrow(pseudo_bulk_seurat_group@meta.data %>% filter(!!as.symbol(comparison) != ident1))
      if (cells1 > 2 & cells2 > 2){
        count_df <- pseudo_bulk_seurat_group@assays$RNA@counts %>% as.data.frame()
        dds <- DESeqDataSetFromMatrix(countData = count_df,
                                      colData = pseudo_bulk_seurat_group@meta.data,
                                      design = formula(paste("~",comparison)))
        dds <- DESeq(dds)
        resDE <- results(dds, name=glue("{comparison}_{ident1}_vs_{ident2}")) %>% as.data.frame() %>% rownames_to_column("ID") %>% filter(!is.na(log2FoldChange)) %>% 
          mutate(!!sym(groups_col) := group)
        
        norm_matrix <- counts(dds, normalized=TRUE)
        
        PathwayActivity_zscore <- progeny(norm_matrix,scale=TRUE, 
                                          organism="Human", top = 100, 
                                          perm = 10000, z_scores = TRUE) %>%
          as.data.frame() %>% rownames_to_column("cell_id") %>% mutate(cell_id = gsub("\\.","-",cell_id)) %>%
          pivot_longer(names_to = "path",values_to = "NES", cols = !cell_id) %>% 
          left_join(pseudo_bulk_seurat_group@meta.data %>% 
                      rownames_to_column("cell_id") %>% 
                      select(cell_id,!!as.symbol(comparison))) %>%
          mutate(path = gsub("\\.","-",path))
        
        library(magrittr)
        pval_pathways <- PathwayActivity_zscore %>% 
          nest(-path) %>%
          mutate(pval = map_dbl(data, ~wilcox.test(
            .x %>% filter(!!as.symbol(comparison) == ident1) %$% NES,
            .x %>% filter(!!as.symbol(comparison) == ident2) %$% NES )$p.value)) %>%
          select(path,pval)
        
        diff_mean_pathways <- rbindlist(lapply(pval_pathways$path, function(path){
          data.frame(path = path,
                     diff_mean = mean(PathwayActivity_zscore %>% 
                                        filter(path == !!path, !!as.symbol(comparison) == ident1) %>% 
                                        pull(NES)) -  mean(PathwayActivity_zscore %>% 
                                                             filter(path == !!path, !!as.symbol(comparison) == ident2) %>% pull(NES)))
        })) %>% left_join(pval_pathways) %>% mutate(!!sym(groups_col) := group)
        
        prog_matrix <- getModel("Human", top = 100) %>% 
          as.data.frame()  %>%
          tibble::rownames_to_column("ID")
        
        pathways_genes <- rbindlist(lapply(unique(diff_mean_pathways$path), function(path){
          weight_matrix <- prog_matrix %>% select(ID,!!as.symbol(path)) %>% filter(!!as.symbol(path) != 0) %>% 
            mutate(w = !!as.symbol(path)) %>% select(-!!as.symbol(path))
          sub_df <- merge(resDE,weight_matrix,by = "ID") %>% mutate(w_score = w * log2FoldChange)
          q <- quantile(sub_df$w_score,0.8)
          sub_df <- sub_df %>% mutate(show_gene = ifelse(w_score >= q,ID,NA),
                                      Pathway = path)
          return(sub_df)
        })) %>% mutate(!!sym(groups_col) := group)
        
        
        TFActivity_zscore <- run_viper(norm_matrix,hs_regulon,
                                       options =  list(minsize = 4, eset.filter = FALSE, 
                                                       cores = 1, verbose = FALSE, nes = TRUE)) %>%
          t() %>%
          as.data.frame() %>% rownames_to_column("cell_id") %>% mutate(cell_id = gsub("\\.","-",cell_id)) %>%
          pivot_longer(names_to = "tf",values_to = "NES", cols = !cell_id) %>% 
          left_join(pseudo_bulk_seurat_group@meta.data %>% 
                      rownames_to_column("cell_id") %>% 
                      select(cell_id,!!as.symbol(comparison))) %>%
          mutate(path = gsub("\\.","_",tf))
        
        library(magrittr)
        pval_tfs <- TFActivity_zscore %>% 
          nest(-tf) %>%
          mutate(pval = map_dbl(data, ~wilcox.test(
            .x %>% filter(!!as.symbol(comparison) == ident1) %$% NES,
            .x %>% filter(!!as.symbol(comparison) == ident2) %$% NES )$p.value)) %>%
          select(tf,pval)
        
        diff_mean_tfs <- rbindlist(lapply(pval_tfs$tf, function(path){
          data.frame(tf = path,
                     diff_mean = mean(TFActivity_zscore %>% 
                                        filter(tf == !!path, !!as.symbol(comparison) == ident1) %>% 
                                        pull(NES)) -  mean(TFActivity_zscore %>% 
                                                             filter(tf == !!path, !!as.symbol(comparison) == ident2) %>% pull(NES)))
        })) %>% left_join(pval_tfs) %>% arrange(pval) %>% mutate(!!sym(groups_col) := group)
        
        tf_genes <- rbindlist(lapply(unique(diff_mean_tfs$tf), function(path){
          weight_matrix <- hs_regulon %>% filter(tf == path) %>% rename("target" = "ID")
          sub_df <- merge(resDE,weight_matrix,by = "ID") %>% mutate(w_score = mor * log2FoldChange)
          q <- quantile(sub_df$w_score,0.8)
          sub_df <- sub_df %>% mutate(show_gene = ifelse(w_score >= q,ID,NA),
                                      TF = path)
          return(sub_df)
        })) %>% mutate(!!sym(groups_col) := group)
        
        
        return(list(resDE,diff_mean_pathways,diff_mean_tfs,pathways_genes,tf_genes))
      } else{
        return(NA)
      }
    }, mc.cores = 10)
    
    diff_expression <- Filter(function(a) any(!is.na(a)), diff_expression)
    
  } else{
    print("Differential expression analysis via meta-analysis")
    ident1 <- pseudo_bulk_seurat@meta.data %>% pull(!!as.symbol(comparison)) %>% unique()
    diff_expression <- lapply(ident1, function(group){
      ident2 <- ident1[ident1 != group]
      diff_expression_ident2 <- mclapply(ident2, function(group2){
        samples <- pseudo_bulk_seurat@meta.data %>% filter(!!as.symbol(comparison) == group | !!as.symbol(comparison) == group2) %>% rownames_to_column("sample_id") %>% pull(sample_id)
        pseudo_bulk_seurat_group <- pseudo_bulk_seurat[,samples]
        cells1 = pseudo_bulk_seurat_group@meta.data %>% filter(!!as.symbol(comparison) == group) %>% rownames_to_column("cell_id")
        cells2 = pseudo_bulk_seurat_group@meta.data %>% filter(!!as.symbol(comparison) == group2) %>% rownames_to_column("cell_id")
        if (nrow(cells1) > 2 & nrow(cells2) > 2){
          count_df <- pseudo_bulk_seurat_group@assays$RNA@counts %>% as.data.frame()
          dds <- DESeqDataSetFromMatrix(countData = count_df,
                                        colData = pseudo_bulk_seurat_group@meta.data,
                                        design = formula(paste("~",comparison)))
          dds <- DESeq(dds)
          logcounts = counts(dds, normalize = T)[,c(cells1$cell_id,cells2$cell_id)]
          objectMA = list(logcounts, c(rep(1,length(cells1$cell_id)),rep(0,length(cells2$cell_id))))
          return(objectMA)
        } else{
          return(NA)
        }
      }, mc.cores = 10)
      
      names(diff_expression_ident2) = ident2
      diff_expression <- Filter(function(a) any(!is.na(a)), diff_expression_ident2)
      library(DExMA)
      res_genes <- metaAnalysisDE(diff_expression, typeMethod = "REM") %>%
        rownames_to_column("ID") %>%
        mutate(!!sym(groups_col) := group)
    })
    diff_expression <- Filter(function(a) any(!is.na(a)), diff_expression)
  }
  
  de <- rbindlist(sapply(diff_expression,"[",1))
  pa <- rbindlist(sapply(diff_expression,"[",2))
  pa_genes <- rbindlist(sapply(diff_expression,"[",4))
  tfa <- rbindlist(sapply(diff_expression,"[",3))
  tfa_genes <- rbindlist(sapply(diff_expression,"[",5))
  return(list("de" = de, 
              "pa" = pa,
              "pa_genes" = pa_genes,
              "tfa" = tfa,
              "tfa_genes" = tfa_genes))
}

pseudo_bulk_seurat <- get_pseudobulk(data_seurat, sample_col = "orig.ident", groups_col = "cluster_name")
all_results <- perform_functional_analysis(pseudo_bulk_seurat_in = pseudo_bulk_seurat,
                                           sample_col = "orig.ident", 
                                           groups_col = "cluster_name", 
                                           comparison = "Condition", 
                                           ident1 = "SLE",
                                           rem = F)

saveRDS(all_results, file = "results/functional_analysis.rds")