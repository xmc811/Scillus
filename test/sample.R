library(Scillus)

# Please install the package first
library(org.Hs.eg.db)
library(fgsea)

hmks <- gmtPathways("~/Downloads/h.all.v7.0.symbols.gmt") # please change the file path

groups <- list.files(path = "./tests/testdata/pdx/") # please change the file path
paths <- paste0("./tests/testdata/pdx/", groups) # please change the file path

scRNA <- purrr::map2(.x = paths, .y = groups, .f = load_scfile, gcol = 2, file = NULL)

plot_qc(scRNA, metrics = "nFeature_RNA")
plot_qc(scRNA, metrics = "nCount_RNA")
plot_qc(scRNA, metrics = "percent.mt")

scRNA <- filter_scdata_list(scRNA, range_mt = c(-Inf, 30), range_nCount = c(500, 50000))

scRNA <- purrr::map(.x = scRNA, .f = NormalizeData)
scRNA <- purrr::map(.x = scRNA, .f = FindVariableFeatures)
scRNA <- purrr::map(.x = scRNA, .f = CellCycleScoring, 
                    s.features = cc.genes$s.genes, 
                    g2m.features = cc.genes$g2m.genes)

plot_qc(scRNA, metrics = "S.Score")
plot_qc(scRNA, metrics = "G2M.Score")

pdx <- IntegrateData(anchorset = FindIntegrationAnchors(object.list = scRNA, dims = 1:30), dims = 1:30)

pdx %<>%
        ScaleData(vars.to.regress = c("nCount_RNA", "percent.mt", "S.Score", "G2M.Score"))
        

pdx %<>%
        RunPCA(npcs = 50, verbose = TRUE)

ElbowPlot(pdx, ndims = 50)

pdx %<>%
        RunUMAP(reduction = "pca", dims = 1:20, n.neighbors = 120) %>%
        FindNeighbors(reduction = "pca", dims = 1:20) %>%
        FindClusters(resolution = 0.2)

plot_scdata(pdx)
plot_scdata(pdx, color_by = "group")
plot_scdata(pdx, split = "group")
plot_scdata(pdx, color_by = "group", split = "group")

plot_stat(pdx, "group_count", tilt_text = TRUE)
plot_stat(pdx, "cluster_count")
plot_stat(pdx, "prop_fill", tilt_text = TRUE)
plot_stat(pdx, "prop_multi", tilt_text = TRUE)

pdx_markers <- FindAllMarkers(pdx, logfc.threshold = 0.1, min.pct = 0, only.pos = T)

plot_heatmap(pdx, pdx_markers, 8)

plot_measure(pdx, features = c("nFeature_RNA", "nCount_RNA", "S.Score", "G2M.Score"), plot_type = "group", meta = TRUE)
plot_measure(pdx, features = c("APOE","MALAT1","TUBB","MMP7","TOE1","AGO4"), plot_type = "cluster")

plot_measure_cluster(dataset = pdx, features = c("APOE","MALAT1","TUBB","MMP7","TOE1","AGO4"))
plot_measure_cluster(dataset = pdx, features = c("APOE","MALAT1","TUBB","MMP7"), plot_type = "group")
plot_measure_cluster(dataset = pdx, features = c("APOE","MALAT1","TUBB","MMP7"), plot_type = "cluster")


plot_all_cluster_go(pdx_markers, ont = "MF")
plot_all_cluster_go(pdx_markers, ont = "CC")
plot_all_cluster_go(pdx_markers, ont = "BP")

plot_GSVA(pdx, hmks)


