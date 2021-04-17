library(Scillus)
library(Seurat)
library(magrittr)
library(ggplot2)
library(tidyverse)

a <- list.files("./test/GSE128531_RAW", full.names = TRUE)
m <- tibble::tibble(file = a, 
                     sample = factor(stringr::str_remove(basename(a), ".csv.gz")),
                     group = factor(rep(c("CTCL", "Normal"), each = 3), levels = c("Normal", "CTCL")))

pal <- tibble::tibble(var = c("sample", "group","seurat_clusters"),
                      pal = c("Set2","Set1","Paired")) 

scRNA <- load_scfile(m)

plot_qc(scRNA, metrics = "percent.mt")
plot_qc(scRNA, metrics = "nFeature_RNA")
plot_qc(scRNA, metrics = "nCount_RNA")
plot_qc(scRNA, metrics = "nCount_RNA", group_by = "group")
plot_qc(scRNA, metrics = "nCount_RNA", group_by = "group", pal_setup = pal)
plot_qc(scRNA, metrics = "nCount_RNA", group_by = "group", pal_setup = c('blue','red'))
plot_qc(scRNA, metrics = "nCount_RNA", pal_setup = c('blue','red','green','yellow','grey','purple'))
plot_qc(scRNA, metrics = "nCount_RNA", pal_setup = c('blue','red','green'))


scRNA <- filter_scdata(scRNA, subset = nFeature_RNA > 500 & percent.mt < 10)

scRNA %<>% 
        purrr::map(.f = NormalizeData) %>%
        purrr::map(.f = FindVariableFeatures) %>%
        purrr::map(.f = CellCycleScoring, 
                   s.features = cc.genes$s.genes, 
                   g2m.features = cc.genes$g2m.genes)

plot_qc(scRNA, metrics = "S.Score", group_by = "group")
plot_qc(scRNA, metrics = "G2M.Score", group_by = "group")

scRNA <- IntegrateData(anchorset = FindIntegrationAnchors(object.list = scRNA, dims = 1:30, k.filter = 50), dims = 1:30)

scRNA %<>%
        ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))

scRNA %<>%
        RunPCA(npcs = 50, verbose = TRUE)

ElbowPlot(scRNA, ndims = 50)

scRNA %<>%
        RunUMAP(reduction = "pca", dims = 1:20, n.neighbors = 30) %>%
        FindNeighbors(reduction = "pca", dims = 1:20) %>%
        FindClusters(resolution = 0.2)

scRNA %<>%
        refactor_seurat(metadata = m)

plot_scdata(scRNA)
plot_scdata(scRNA, pal_setup = pal)
plot_scdata(scRNA, pal_setup = "Set3")
plot_scdata(scRNA, 
            pal_setup = c('red','yellow','blue','green','cyan','purple','orange'))
plot_scdata(scRNA, color_by = "sample", 
            pal_setup = c('red','yellow','blue','green','cyan','purple','orange'))
plot_scdata(scRNA, split_by = "sample")
plot_scdata(scRNA, split_by = "No Split")
plot_scdata(scRNA, color_by = "group", split_by = "seurat_clusters", pal_setup = pal)

plot_stat(scRNA, "group_count", group_by = "seurat_clusters", pal_setup = pal)
plot_stat(scRNA, "group_count", group_by = "sample")
plot_stat(scRNA, "group_count", group_by = "group", pal_setup = c('red','blue'))

plot_stat(scRNA, "prop_fill", group_by = "group", 
          pal_setup = c('red','yellow','blue','green','cyan','purple','orange'))
plot_stat(scRNA, "prop_multi", group_by = "sample", pal_setup = pal)

markers <- FindAllMarkers(scRNA, logfc.threshold = 0.1, min.pct = 0, only.pos = T)

plot_heatmap(dataset = scRNA, 
              markers = markers,
              sort_var = c("seurat_clusters","sample"),
              anno_var = c("seurat_clusters","sample","percent.mt","S.Score","G2M.Score"),
              anno_colors = list("Set2",
                                 c('red','orange','yellow','purple','blue','green'),
                                 "Reds",
                                 c('blue','white','red'),
                                 "Greens"),
             hm_limit = c(-1.5,0,1.5),
             hm_colors = c('purple','black','yellow'))


plot_cluster_go(markers, cluster_name = '1', org = "human", ont = "CC")

plot_all_cluster_go(markers, org = 'human', ont = "CC")

plot_measure(dataset = scRNA, measures = c("KRT14","S100A8","FAM138A","percent.mt"), 
             group_by = "group", 
             pal_setup = pal)

plot_measure_dim(dataset = scRNA, 
                 measures = c("nFeature_RNA","nCount_RNA","percent.mt","KRT14"))

plot_measure_dim(dataset = scRNA, 
                 measures = c("nFeature_RNA","nCount_RNA","percent.mt","KRT14"),
                 split_by = "No Split")

plan("multiprocess", workers = 2)
de <- find_diff_genes(dataset = scRNA_1, 
                      clusters = as.character(0:6),
                      comparison = c("group", "CTCL", "Normal"),
                      logfc = 0)

gsea_res <- test_GSEA(de, 
                      pathway = pathways.hallmark)

plot_GSEA(gsea_res)


# test

ha = HeatmapAnnotation(a = 1:10)
names(ha) <- 'foo'
ha@anno_list[[1]]@color_mapping@name <- 'foo'
ha@anno_list[[1]]@name_param$label <- 'foo'
names(ha@anno_list)

