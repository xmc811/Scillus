library(magrittr)
library(Seurat)

a1 <- list.files("./test/GSE128531_RAW", full.names = TRUE)
m1 <- tibble::tibble(file = a1, sample = stringr::str_remove(basename(a1), ".csv.gz"))

scRNA_1 <- load_scfile(m1)

plot_qc(scRNA_1, metrics = "percent.mt")
plot_qc(scRNA_1, metrics = "nFeature_RNA")
plot_qc(scRNA_1, metrics = "nCount_RNA")


scRNA_1 <- filter_scdata(scRNA_1, subset = nFeature_RNA > 500 & percent.mt < 10)

scRNA_1 %<>% 
        purrr::map(.f = NormalizeData) %>%
        purrr::map(.f = FindVariableFeatures) %>%
        purrr::map(.f = CellCycleScoring, 
                   s.features = cc.genes$s.genes, 
                   g2m.features = cc.genes$g2m.genes)

plot_qc(scRNA_1, metrics = "S.Score")
plot_qc(scRNA_1, metrics = "G2M.Score")

scRNA_1 <- IntegrateData(anchorset = FindIntegrationAnchors(object.list = scRNA_1, dims = 1:30, k.filter = 50), dims = 1:30)

scRNA_1 %<>%
        ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))

scRNA_1 %<>%
        RunPCA(npcs = 50, verbose = TRUE)

ElbowPlot(scRNA_1, ndims = 50)

scRNA_1 %<>%
        RunUMAP(reduction = "pca", dims = 1:20, n.neighbors = 30) %>%
        FindNeighbors(reduction = "pca", dims = 1:20) %>%
        FindClusters(resolution = 0.2)

plot_scdata(scRNA_1)
plot_scdata(scRNA_1, color_by = "sample")
plot_scdata(scRNA_1, split = "sample")
plot_scdata(scRNA_1, color_by = "sample", split = "sample")

plot_stat(scRNA_1, "sample_count")
plot_stat(scRNA_1, "cluster_count")
plot_stat(scRNA_1, "prop_fill")
plot_stat(scRNA_1, "prop_multi")

markers_1 <- FindAllMarkers(scRNA_1, logfc.threshold = 0.1, min.pct = 0, only.pos = T)

plot_heatmap(dataset = scRNA_1, 
              markers = markers_1,
              sort_var = c("seurat_clusters","sample"),
              anno_var = c("seurat_clusters", "sample","percent.mt","S.Score","G2M.Score"),
              anno_colors = c("Set2","Dark2","Reds","Blues","Greens"))



