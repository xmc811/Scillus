
# Integrative analysis
DefaultAssay(object = gfp_combined) <- "integrated"

gfp_combined <- subset(gfp_combined, cells = setdiff(names(gfp_combined@active.ident), c(doublets[[1]],doublets[[3]])))

gfp_combined <- ScaleData(object = gfp_combined,
                          vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"),
                          verbose = T)

gfp_combined <- analyze_merged(gfp_combined, group.levels = stages,
                               verbose = T, npcs = 50, dims = 1:20, nnei = 70, k.param = 70, min.dist = 0.8, spread = 1.5, resolution = 0.2)

# Visualization

svglite(file = "./figures/gfp_merge.svg")
plot_merge(gfp_combined)
dev.off()

plot_cluster(gfp_combined, label = F)

plot_split(gfp_combined, colors = get_colors(seq(levels(gfp_combined$seurat_clusters)), pal = "Set3"))

svglite(file = "./figures/gfp_markers.svg", width = 15, height = 8)
plot_features(gfp_combined, features = c("Krt14","Csn3","Lyz2","Ptn","S100a9","Col1a1","Cd3d"), ncol = 4)
dev.off()

# Identify cell markers

gfp_markers <- FindAllMarkers(gfp_combined, only.pos = T, logfc.threshold = 0.1)

# Heatmap

plot_heatmap(gfp_combined, gfp_markers, 6, cluster_pal = c("Set3"))

# Relabeling

gfp_labels <- c("Basal 1", 
                "Luminal 1", 
                "Macrophage", 
                "Luminal 2",
                "Neutrophil",
                "NS Immune",
                "Basal 2", 
                "Fibroblast")

gfp_levels <- gfp_labels[c(1,7,2,4,3,5,6,8)]

gfp_combined <- rename_cluster(gfp_combined, gfp_labels)

# New visualization

gfp_index <- c(5:8,4,2,10,12)

svglite(file = "./figures/gfp_cluster.svg")
plot_cluster(gfp_combined, label = F, levels = gfp_levels, self_set_color = T, self_colors = get_colors(gfp_index))
dev.off()

# Statistics

svglite(file = "./figures/gfp_group_count.svg")
plot_stat(gfp_combined, "group_count", group_levels = stages, cluster_levels = gfp_levels, plot_ratio = 3)
dev.off()

svglite(file = "./figures/gfp_cluster_count.svg")
plot_stat(gfp_combined, "cluster_count", group_levels = stages, cluster_levels = gfp_levels,
          self_set_color = T, self_colors = get_colors(gfp_index))
dev.off()

svglite(file = "./figures/gfp_prop_fill.svg")
plot_stat(gfp_combined, "prop_fill", group_levels = stages, cluster_levels = gfp_levels, plot_ratio = 3,
          self_set_color = T, self_colors = get_colors(gfp_index))
dev.off()

svglite(file = "./figures/gfp_prop_diverge.svg")
plot_stat(gfp_combined, "prop_diverge", group_levels = stages, cluster_levels = gfp_levels, plot_ratio = 0.8)
dev.off()

plot_stat(gfp_combined, "prop_multi", group_levels = stages, cluster_levels = gfp_levels, plot_ratio = 3)


# DE analysis

gfp_diff <- find_diff_genes(dataset = gfp_combined, clusters = gfp_levels, groups = stages, logfc = 0)

gfp_hallmark_gsea <- test_GSEA(gfp_diff, gfp_levels, pathways.hallmark)

svglite(file = "./figures/gfp_gsea.svg", width = 12, height = 12)
plot_GSEA(gfp_hallmark_gsea, p_cutoff = 0.1, levels = gfp_levels)
dev.off()


# Test

head(Idents(gfp_combined))



test <- FindMarkers(gfp_combined, 
            ident.1 = "Basal 1",
            ident.2 = "Basal 2",
            logfc.threshold = 0,
            assay = "RNA")
test %>%
            add_column(feature = rownames(test), .before = 1) %>%
            add_column(cluster = "Basal 1 vs 2", .after = 1) %>%    
            as_tibble() %>%
            test_GSEA(clusters = c("Basal 1 vs 2"), pathway = pathways.hallmark) %>%
            plot_GSEA(p_cutoff = 0.1, levels = c("Basal 1 vs 2"))


cluster_diff <- function(object, clusters, thr) {
        
        test <- FindMarkers(object, 
                            ident.1 = clusters[1],
                            ident.2 = clusters[2],
                            logfc.threshold = thr,
                            assay = "RNA")
        
        label <- paste0(clusters[1], clusters[2])
        
        test %>%
                add_column(feature = rownames(test), .before = 1) %>%
                add_column(cluster = label, .after = 1) %>%    
                as_tibble() %>%
                test_GSEA(clusters = label, pathway = pathways.hallmark) %>%
                plot_GSEA(p_cutoff = 0.05, levels = label)
        
}

cluster_diff(gfp_combined, c("Basal 1", "Basal 2"), thr = 0)

cluster_diff(gfp_combined, c("Basal 1", "Luminal 1"), thr = 0)

cluster_diff(gfp_combined, c("Basal 2", "Luminal 1"), thr = 0)

cluster_diff(gfp_combined, c("Luminal 2", "Luminal 1"), thr = 0)


test <- gfp_combined

test$celltype.group <- paste(Idents(object = test), test$group, sep = "_")
Idents(object = test) <- "celltype.group"

head(Idents(test))

cluster_diff(test, c("Basal 1_Early", "Basal 2_Early"), thr = 0)
cluster_diff(test, c("Basal 1_Advanced", "Basal 2_Advanced"), thr = 0)

cluster_diff(test, c("Basal 1_Early", "Luminal 1_Early"), thr = 0)
cluster_diff(test, c("Basal 1_Advanced", "Luminal 1_Advanced"), thr = 0)

cluster_diff(test, c("Basal 2_Early", "Luminal 1_Early"), thr = 0)
cluster_diff(test, c("Basal 2_Advanced", "Luminal 1_Advanced"), thr = 0)

plot_GSEA(test_GSEA(gfp_diff, gfp_levels, pathways.hallmark[1]), levels = gfp_levels)

pathways.IL17 <- list()
pathways.IL17$HALLMARK_IL17 <- pathways.hallmark$HALLMARK_TNFA_SIGNALING_VIA_NFKB

plot_GSEA(test_GSEA(gfp_diff, gfp_levels, pathways.IL17), levels = gfp_levels)
plot_GSEA(test_GSEA(cd45_diff, cd45_levels, pathways.IL17), levels = cd45_levels)

tibble(x = cd45_levels, y = rnorm(n = 11, mean = 2, sd = 1), z = runif(n = 11, min = 1, max = 2.5), p = 1) %>%
        ggplot(aes(x = factor(p), y = factor(x))) + 
        geom_point(aes(size = abs(z), color = y)) +
        scale_size(name = "Normalized\nEnrichment\nScore Size") +
        scale_color_gradient2(name = bquote(-log[10]~"Adj. p-value"), low = 'dodgerblue1', mid = 'grey', high = 'red', midpoint = 0) +
        coord_flip() +
        geom_point(aes(shape = factor(p)), size = 5.5, stroke = 1) +
        scale_shape_manual(name = "Adj. p-value", values=c(NA, 0)) +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank())


test_GSEA(gfp_diff, gfp_levels, pathway = pathways.IL17)



