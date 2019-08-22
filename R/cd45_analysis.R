
# Integrative analysis
DefaultAssay(object = cd45_combined) <- "integrated"

cd45_combined <- subset(cd45_combined, cells = setdiff(names(cd45_combined@active.ident), c(doublets[[2]],doublets[[4]])))

cd45_combined <- ScaleData(object = cd45_combined,
                          vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"),
                          verbose = T)

cd45_combined <- analyze_merged(cd45_combined, group.levels = stages, 
                                verbose = T, npcs = 50, dims = 1:30, nnei = 70, k.param = 70, min.dist = 0.6, spread = 1.5, resolution = 0.3)

# Visualization

svglite(file = "./figures/cd45_merge.svg")
plot_merge(cd45_combined)
dev.off()

plot_cluster(cd45_combined, label = F)

plot_split(cd45_combined, colors = get_colors(seq(levels(cd45_combined$seurat_clusters)), pal = "Set3"))

svglite(file = "./figures/cd45_markers.svg", width = 15, height = 8)
plot_features(cd45_combined, features = c("Cd3d","Gzma","S100a8","Ebf1","Lyz2","Siglech","Iglv1","Col1a1","Plbd1","Ptprc"), ncol = 5)
dev.off()

svglite(file = "./figures/cd45_markers_2.svg", width = 15, height = 8)
plot_features(cd45_combined, features = c("Cd4","Cd8a", "Cd8b1"), ncol = 3)
dev.off()

# Identify cell markers

cd45_markers <- FindAllMarkers(cd45_combined, only.pos = T, logfc.threshold = 0.1)

# Heatmap

plot_heatmap(cd45_combined, cd45_markers, 5, cluster_pal = c("Set3"))

# Relabeling

cd45_labels <- c("T Cell 1",
                 "T Cell 2",
                 "NK Cell", 
                 "Neutrophil", 
                 "B Cell",
                 "Macrophage",
                 "pDC",
                 "Plasma Cell",
                 "cDC",
                 "T Cell 3",
                 "T Cell 4"
                 )

cd45_levels <- cd45_labels[c(1,2,10,11,3,5,8,6,4,7,9)]

cd45_combined <- rename_cluster(cd45_combined, cd45_labels)

# New visualization

cd45_index <- c(1:4,12,10,9,8,6,7,5)

svglite(file = "./figures/cd45_cluster.svg")
plot_cluster(cd45_combined, label = F, levels = cd45_levels, self_set_color = T, self_colors = get_colors(cd45_index))
dev.off()

# Statistics

svglite(file = "./figures/cd45_group_count.svg")
plot_stat(cd45_combined, "group_count", group_levels = stages, cluster_levels = cd45_levels, plot_ratio = 3)
dev.off()

svglite(file = "./figures/cd45_cluster_count.svg")
plot_stat(cd45_combined, "cluster_count", group_levels = stages, cluster_levels = cd45_levels,
          self_set_color = T, self_colors = get_colors(cd45_index))
dev.off()

svglite(file = "./figures/cd45_prop_fill.svg")
plot_stat(cd45_combined, "prop_fill", group_levels = stages, cluster_levels = cd45_levels, plot_ratio = 3,
          self_set_color = T, self_colors = get_colors(cd45_index))
dev.off()

svglite(file = "./figures/cd45_prop_diverge.svg")
plot_stat(cd45_combined, "prop_diverge", group_levels = stages, cluster_levels = cd45_levels, plot_ratio = 0.8)
dev.off()

plot_stat(cd45_combined, "prop_multi", group_levels = stages, cluster_levels = cd45_levels, plot_ratio = 1.5)


# DE analysis

cd45_diff <- find_diff_genes(dataset = cd45_combined, clusters = cd45_levels, groups = stages, logfc = 0)

cd45_hallmark_gsea <- test_GSEA(cd45_diff, cd45_levels, pathways.hallmark)

svglite(file = "./figures/cd45_gsea.svg", width = 12, height = 12)
plot_GSEA(cd45_hallmark_gsea, p_cutoff = 0.1, levels = cd45_levels)
dev.off()



# Test

test2 <- cd45_combined

test2$celltype <- Idents(test2)

DefaultAssay(test2) <- "RNA"

svglite(file = "./figures/IL17_1.svg", width = 12, height = 12)
plots <- VlnPlot(test2, features = c("S100a9", "S100a8", "Fosb"), split.by = "group", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()

svglite(file = "./figures/IL17_2.svg", width = 12, height = 12)
plots <- VlnPlot(test2, features = c("Il17a", "Il17ra"), split.by = "group", group.by = "celltype", 
                 pt.size = 0, combine = FALSE)
CombinePlots(plots = plots, ncol = 1)
dev.off()




