
# Neutrophil subset

nf_combined <- subset(x = cd45_combined, idents = c("Neutrophil"))

nf_combined <- ScaleData(object = nf_combined,
                            vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"),
                            verbose = T)

nf_combined <- analyze_merged(nf_combined, group.levels = stages,
                                 verbose = T, npcs = 50, dims = 1:20, nnei = 25, k.param = 25, min.dist = 0.05, spread = 0.5, resolution = 0.1)

nf_combined <- subset(x = nf_combined, idents = as.character(0:4))

nf_combined <- analyze_merged(nf_combined, group.levels = stages,
                                 verbose = T, npcs = 50, dims = 1:20, nnei = 40, k.param = 40, min.dist = 0.005, spread = 0.5, resolution = 0.2)

# Visualization

svglite(file = "./figures/nf_merge.svg")
plot_merge(nf_combined)
dev.off()

svglite(file = "./figures/nf_cluster.svg")
plot_cluster(nf_combined, label = F)
dev.off()

plot_split(nf_combined, colors = get_colors(seq(levels(nf_combined$seurat_clusters)), pal = "Set3"))

svglite(file = "./figures/nf_markers.svg", width = 15, height = 8)
plot_features(nf_combined, features = c("Cd4","Cd8a","Cd8b1","Satb1","Ctla4","Pdcd1","Cd3d"), ncol = 4)
dev.off()

# Identify cell markers

nf_markers <- FindAllMarkers(nf_combined, only.pos = T, logfc.threshold = 0.1)

# Heatmap

plot_heatmap(nf_combined, nf_markers, 10, cluster_pal = c("Set3"))

# Relabeling

tcell_labels <- c("CD8+ Non-Ex",
                  "CD8+ Non-Ex",
                  "CD4+ Non-Ex",
                  "CD4+ Ex",
                  "CD8+ Non-Ex",
                  "CD8+ Ex",
                  "CD4+ Ex",
                  "CD8+ Ex",
                  "CD8+ Ex",
                  "CD4+ Ex")

tcell_levels <- tcell_labels[c(3,4,1,6)]

tcell_combined <- rename_cluster(tcell_combined, tcell_labels)

tcell_index <- 1:4

# Statistics

svglite(file = "./figures/tcell_group_count.svg")
plot_stat(tcell_combined, "group_count", group_levels = stages, cluster_levels = tcell_levels, plot_ratio = 3)
dev.off()

svglite(file = "./figures/tcell_cluster_count.svg")
plot_stat(tcell_combined, "cluster_count", group_levels = stages, cluster_levels = tcell_levels,
          self_set_color = T, self_colors = get_colors(tcell_index))
dev.off()

svglite(file = "./figures/tcell_prop_fill.svg")
plot_stat(tcell_combined, "prop_fill", group_levels = stages, cluster_levels = tcell_levels, plot_ratio = 3,
          self_set_color = T, self_colors = get_colors(tcell_index))
dev.off()

svglite(file = "./figures/tcell_prop_diverge.svg")
plot_stat(tcell_combined, "prop_diverge", group_levels = stages, cluster_levels = tcell_levels, plot_ratio = 0.8)
dev.off()

plot_stat(tcell_combined, "prop_multi", group_levels = stages, cluster_levels = tcell_levels, plot_ratio = 3)

# T cell exhaustion score

exhaust_genes <- read.table("./refs/exhaustion_gene_hs", stringsAsFactors = F)[[1]]

tcell_combined <- add_program_score(tcell_combined, features = exhaust_genes, org = "mouse", nbin = 20, ctrl = 10, name = "Exhaust_Score")

svglite(file = "./figures/tcell_ex_score.svg")
plot_measure(tcell_combined, measure = "Exhaust_Score1", plot_type = "cluster", group_levels = stages, cluster_levels = tcell_levels)
dev.off()

plot_measure(tcell_combined, measure = "Exhaust_Score1", plot_type = "cluster_group", group_levels = stages, cluster_levels = tcell_levels)


# Test