
# Merging data

wahl_combined <- IntegrateData(anchorset = FindIntegrationAnchors(object.list = c(scRNA[c(1,3)], scRNA_wahl[1:4]), dims = 1:30), dims = 1:30)


# Integrative analysis
DefaultAssay(object = wahl_combined) <- "integrated"

wahl_combined <- ScaleData(object = wahl_combined,
                          vars.to.regress = c("percent.mt","nCount_RNA"),
                          verbose = T)

wahl_combined <- analyze_merged(wahl_combined, group.levels = c(stages, wahl_stages),
                               verbose = T, npcs = 50, dims = 1:50, nnei = 110, k.param = 110, min.dist = 0.8, spread = 1.5, resolution = 0.5)


# Subset

epi <- subset(wahl_combined, subset = seurat_clusters %in% c(0:5, 7:8))
plot_split(epi, colors = get_colors(c(6,5,12,2,8,1,10,9), pal = "Paired"))

# Visualization

plot_merge(wahl_combined, colors = get_colors(seq(levels(wahl_combined$group)), "Set1"))

plot_cluster(wahl_combined, label = F, self_set_color = T, self_colors = get_colors(seq(levels(wahl_combined$seurat_clusters)), pal = "Paired"))

plot_split(wahl_combined, colors = get_colors(seq(levels(wahl_combined$seurat_clusters)), pal = "Paired"))

plot_features(wahl_combined, features = c("Krt14","Csn3","Lyz2","Ptn","S100a9","Col1a1","Cd3d","Krt5"), ncol = 4)


# Density plot

plot_density <- function(data, stage, scale_ind) {
        
        cells.use <- colnames(data)[which(data[[]]['group'] == stage)]
        
        obj <- subset(data, cells = cells.use)
        coor <- as_tibble(obj@reductions$umap@cell.embeddings)
       
        ggplot(coor, aes(x = UMAP_1, y = UMAP_2)) +
                stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
                scale_fill_gradient2(low = "white", mid = scales[[scale_ind]][1], high = scales[[scale_ind]][2]) +
                scale_x_continuous(limits = c(-12,2), expand = c(0, 0)) +
                scale_y_continuous(limits = c(-20,20), expand = c(0, 0)) +
                ggtitle(stage) +
                theme(panel.border = element_rect(colour = "black", fill=NA, size=1, linetype = 1),
                      plot.title = element_text(face = 'plain'),
                      legend.position = "none",
                      axis.line=element_blank(),
                      axis.title.x=element_blank(),
                      axis.title.y=element_blank(),
                      aspect.ratio = 1)
        
}

scales <- list(c("#fff7ec", "#b30000"),
               c("#f7fbff", "#08306b"),
               c("#f7fcf5", "#00441b"),
               c("#f7fcfd", "#4d004b")
               )

p_den <- purrr::map2(.x = c(stages, wahl_stages), .y = 4, .f = plot_density, data = epi)

plot_grid(plotlist = p_den, ncol = 2)

purrr::map(.x = p_den, plot)

# Test


