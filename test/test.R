library(Scillus)

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
plot_qc(scRNA, metrics = "nCount_RNA", group_by = "group", pal_setup = pal)


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

plot_scdata(scRNA, pal_setup = pal)
plot_scdata(scRNA, color_by = "sample")
plot_scdata(scRNA, split_by = "sample")
plot_scdata(scRNA, color_by = "group", split_by = "seurat_clusters", pal_setup = pal)

plot_stat(scRNA, "group_count", group_by = "group")
plot_stat(scRNA, "cluster_count", group_by = "group")
plot_stat(scRNA, "prop_fill", group_by = "group")
plot_stat(scRNA, "prop_multi", group_by = "group")

markers_1 <- FindAllMarkers(scRNA_1, logfc.threshold = 0.1, min.pct = 0, only.pos = T)

plot_heatmap(dataset = scRNA_1, 
              markers = markers_1,
              sort_var = c("seurat_clusters","sample"),
              anno_var = c("seurat_clusters", "sample","percent.mt","S.Score","G2M.Score"),
              anno_colors = c("Set2","Dark2","Reds","Blues","Greens"))


plot_cluster_go(markers_1, cluster_name = '1', org = "human", ont = "CC")

plot_all_cluster_go(markers_1, org = 'human', ont = "CC")

plot_measure(dataset = scRNA_1, measures = c("KRT14","S100A8","FAM138A","percent.mt"), 
             group_by = "group", 
             pal_setup = pal)

plot_measure_dim(dataset = scRNA_1, measures = c("nFeature_RNA","nCount_RNA","percent.mt"))

plan("multiprocess", workers = 2)
de <- find_diff_genes(dataset = scRNA_1, 
                      clusters = as.character(0:6),
                      comparison = c("group", "CTCL", "Normal"),
                      logfc = 0)

gsea_res <- test_GSEA(de, 
                      pathway = pathways.hallmark)

plot_GSEA(gsea_res)
