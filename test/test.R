
library(Seurat)
library(tidyverse)
library(data.table)

a1 <- list.files("./test/GSE128531_RAW", full.names = TRUE)

m1 <- tibble(file = a1,
            sample = str_remove(basename(a1), ".csv.gz"))

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

plot_heatmap(scRNA_1, markers_1, 8)

pdf("./sample/heatmap.pdf", width = 17, height = 8)
plot_heatmap(scRNA_2, markers_2, 5)
dev.off()



genes <- get_top_genes(scRNA_1, markers_1, 8)

genes <- markers_1 %>%
        arrange(desc(avg_logFC)) %>%
        head(100) %>%
        `[[`('gene')

library(ComplexHeatmap)


plot_heatmap <- function(dataset, 
                          markers,
                          sort_var = c('seurat_clusters', 'sample'),
                          n = 8, 
                          anno_var, 
                          anno_colors) {
        
        mat <- GetAssayData(object = dataset, assay = "integrated", slot = "scale.data")
        genes <- get_top_genes(dataset, markers, n)
        
        mat <- mat[match(genes, rownames(mat)),]
        
        anno <- dataset@meta.data %>%
                rownames_to_column(var = "barcode") %>%
                arrange(!!!syms(sort_var))
        
        mat <- t(mat)
        mat <- mat[match(anno$barcode, rownames(mat)),]
        mat <- t(mat)
        
        
        annos <- list()
        
        for (i in seq_along(1:length(anno_var))) {
                
                value <- anno[[anno_var[i]]]
                
                if (is.numeric(value)) {
                        
                        n <- brewer.pal.info[anno_colors[i],]['maxcolors'][[1]]
                        pal <- brewer.pal(n = n, name = anno_colors[i])
                        
                        col_fun <- colorRamp2(c(min(value), median(value), max(value)), 
                                              c(pal[2], pal[(n+1)/2], pal[n-1]))
                        
                        ha <- HeatmapAnnotation(a = anno[[anno_var[i]]],
                                                col = list(a = col_fun),
                                                border = TRUE)
                } else {
                        
                        l <- levels(factor(anno[[anno_var[i]]]))
                        col <- get_palette(ncolor = length(l), 
                                           palette = anno_colors[i])
                        names(col) <- l
                        col <- list(a = col)
                        
                        ha <- HeatmapAnnotation(a = anno[[anno_var[i]]],
                                                col = col,
                                                border = TRUE)
                }
                
                names(ha) <- anno_var[i]
                names(ha@anno_list) <- anno_var[i]
                ha@anno_list[[1]]@color_mapping@name <- anno_var[i]
                ha@anno_list[[1]]@name_param$label <- anno_var[i]
                
                annos[[i]] <- ha
        }
        
        annos <- do.call(c, annos)
        
        annos@gap <- rep(unit(1,"mm"), length(annos))
        
        ht <- Heatmap(mat,
                      cluster_rows = FALSE,
                      cluster_columns = FALSE,
                      heatmap_legend_param = list(direction = "horizontal",
                                                  legend_width = unit(6, "cm"),
                                                  title = "Expression"),
                      col = colorRamp2(c(-2, 0, 2), c("#4575b4", "white", "#d73027")),
                      show_column_names = FALSE,
                      row_names_side = "left",
                      top_annotation = annos)
        
        draw(ht, 
             heatmap_legend_side = "bottom",
             annotation_legend_side = "right")
        
}

plot_heatmap2(dataset = scRNA_1, 
              markers = markers_1,
              sort_var = c("seurat_clusters","sample"),
              anno_var = c("seurat_clusters", "sample","percent.mt","S.Score","G2M.Score"),
              anno_colors = c("Set2","Dark2","Reds","Blues","Greens"))

library(circlize)


ha@gap <- unit(2,"mm")

length(ha@gap)


