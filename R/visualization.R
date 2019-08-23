
# These functions are used for visualization

plot_qc <- function(data_list, metrics) {
        
        qc <- list()
        
        for (i in seq(length(data_list))) {
                qc[[i]] <- tibble(value = data_list[[i]][[metrics]][[1]],
                                  sample = names(data_list)[i])
        }
        qc <- do.call(rbind, qc)
        ggplot(qc) + 
                geom_boxplot(aes(x = sample, y = value, fill = sample)) +
                scale_fill_brewer(palette = "Set3") +
                labs(y = metrics)
}

plot_merge <- function(dataset, reduction = "umap", group.by = "group", 
                       colors = c('#92c5de','#d6604d'), legend.title = "Group", labels = levels(dataset$group)) {
        
        p <- DimPlot(object = dataset, reduction = reduction, group.by = group.by)
        
        p  + scale_color_manual(values = colors, name = legend.title, labels = labels) +
                theme(panel.border = element_rect(colour = "black", fill=NA, size=1, linetype = 1),
                      axis.line=element_blank(),
                      aspect.ratio = 1
                )
}

plot_split <- function(dataset, reduction = "umap", split.by = "group", 
                       colors = c('#92c5de','#d6604d'), legend.title = "Cluster", labels = levels(dataset$seurat_clusters)) {
        
        p <- DimPlot(object = dataset, reduction = reduction, split.by = split.by)
        
        p  + scale_color_manual(values = colors, name = legend.title, labels = labels) +
                theme(panel.border = element_rect(colour = "black", fill = NA, size = 1, linetype = 1),
                      strip.text.x = element_text(face = "plain", vjust = 1),
                      axis.line=element_blank(),
                      aspect.ratio = 1
                )
}

plot_cluster <- function(dataset, reduction = "umap", label = T, levels = NULL,
                         self_set_color = F,
                         self_colors,
                         palette = c("Set2", "Paired")) {
        
        ncolor <- length(levels(Idents(dataset)))
        colors <- if (self_set_color) self_colors else (get_palette(ncolor, palette))
        
        tmp <- dataset
        
        if (is.null(levels) == F) {
                Idents(tmp) <- factor(Idents(tmp), levels = levels)
        }
        
        p <- DimPlot(object = tmp, reduction = reduction, label = label)
        p + scale_color_manual(values = colors, name = "Clusters") +
                theme(panel.border = element_rect(colour = "black", fill=NA, size=1, linetype = 1),
                      axis.line=element_blank(),
                      strip.text = element_blank(),
                      aspect.ratio = 1
                )
}

plot_features <- function(dataset, features, ncol) {
        
        DefaultAssay(dataset) <- "RNA"
        p_gene <- FeaturePlot(object = dataset, 
                              features = features, 
                              min.cutoff = "q9",
                              cols = rev(colormap(colormap = colormaps$density, nshades = 72, format = "hex",
                                                  alpha = 1, reverse = FALSE)), combine = F)
        
        p_gene <- lapply(X = p_gene, 
                         FUN = function(x) 
                                 x + theme(plot.title = element_text(face = 'plain'),
                                           panel.border = element_rect(colour = "black", fill=NA, size=1, linetype = 1),
                                           axis.line=element_blank(),
                                           axis.title.x=element_blank(),
                                           axis.title.y=element_blank(),
                                           legend.position = 'none',
                                           aspect.ratio = 1)
        )
        CombinePlots(plots = p_gene, ncol = ncol)
}

plot_heatmap <- function(dataset, markers, nfeatures,
                         cluster_pal = c("Paired", "Set2", "Set1"),
                         group_colors = c('#92c5de','#d6604d')) {
        
        df <- as_tibble(cbind(colnames(dataset), dataset$seurat_clusters, dataset$group))
        colnames(df) <- c("barcode","cluster","group")
        df$cluster <- as.numeric(df$cluster)
        df %<>%
                arrange(cluster, group)
        
        genes <- get_top_genes(dataset, markers, nfeatures)
        
        p_heat <- DoHeatmap(object = dataset, assay = 'integrated', features = genes, 
                            group.bar = F, cells = df$barcode, raster = F, draw.lines = F)
        
        p_pos_y <- ggplot_build(plot = p_heat)$layout$panel_params[[1]]$y.range
        
        ncol <- length(levels(Idents(dataset)))
        
        pal1 <- get_palette(ncolor = ncol, palette = cluster_pal)
        col1 <- pal1[as.numeric(df$cluster)]
        
        pal2 <- group_colors
        col2 <- pal2[as.numeric(factor(df$group))]
        
        p_heat + 
                annotation_raster(t(col2), -Inf, Inf, max(p_pos_y)+0.5, max(p_pos_y)+1.5) +
                annotation_raster(t(col1), -Inf, Inf, max(p_pos_y)+2, max(p_pos_y)+3) +
                coord_cartesian(ylim = c(0, max(p_pos_y)+4), clip = 'off') +
                scale_fill_gradient2(low = '#377eb8', high = '#e41a1c', mid = 'white', midpoint = 0) +
                guides(colour="none")
        
}

plot_stat <- function(dataset, plot_type, 
                      group_levels, cluster_levels,
                      self_set_color = F,
                      self_colors,
                      group_colors = c('#92c5de','#d6604d'),
                      palette = c("Set3", "Paired"),
                      plot_ratio = 1,
                      text_size = 10) {
        
        stat <- as_tibble(cbind(group = as.character(dataset$group), cluster = as.character(Idents(dataset))))
        stat %<>%
                mutate(group = factor(group, levels = group_levels),
                       cluster = factor(cluster, levels = cluster_levels)) %>%
                group_by(group, cluster) %>%
                summarise(n = n()) %>%
                mutate(freq = n / sum(n))
        
        cluster_colors <- if(self_set_color) self_colors else (get_palette(length(levels(Idents(dataset))), palette = palette))
        
        thm <- theme(aspect.ratio = plot_ratio,
                     legend.title = element_text(size = text_size),
                     legend.text = element_text(size = text_size),
                     axis.title = element_text(size = text_size),
                     axis.text = element_text(size = text_size),
                     axis.title.x = element_blank()
        )
        
        switch(plot_type,
               group_count = stat %>%
                       group_by(group) %>%
                       summarise(sum(n)) %>%
                       ggplot() +
                       geom_col(aes(x = group, y = `sum(n)`, fill = group)) +
                       geom_text(aes(x = group, y = `sum(n)`, label = `sum(n)`), 
                                 vjust = -0.5, size = text_size * 0.35) +
                       scale_fill_manual(values = group_colors, name = "Group") + 
                       labs(y = "Counts") + thm,
               
               cluster_count = stat %>%
                       group_by(cluster) %>%
                       summarise(sum(n)) %>%
                       ggplot() +
                       geom_col(aes(x = cluster, y = `sum(n)`, fill = cluster)) +
                       geom_text(aes(x = cluster, y = `sum(n)`, label = `sum(n)`), 
                                 vjust = -0.5, size = text_size * 0.35) +
                       scale_fill_manual(values = cluster_colors, name = "Cluster") + 
                       labs(y = "Counts") + 
                       theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
                       thm,
               
               prop_fill = ggplot(stat) + 
                       geom_bar(aes(x = group, y = freq, fill = cluster), position = "fill", stat = "identity") +
                       scale_y_continuous(labels = scales::percent) +
                       scale_fill_manual(values = cluster_colors, name = "Cluster") +
                       labs(y = "Proportion") + thm,
               
               prop_diverge = stat %>%
                       mutate(cluster = fct_rev(cluster)) %>%
                       mutate(freq = ifelse(group == group_levels[1], -freq, freq)) %>%
                       mutate(freq = round(freq, 3)) %>%
                       ggplot() + 
                       geom_bar(aes(x=cluster, y = freq, fill = group), stat = "identity") +
                       geom_text(aes(x=cluster, y = freq + 0.03 * sign(freq), label = percent(abs(freq), digits = 1)), size = text_size * 0.35) +
                       coord_flip() +
                       scale_fill_manual(values = group_colors, name = "Group") +
                       scale_y_continuous(breaks = pretty(c(stat$freq, -stat$freq)),
                                          labels = scales::percent(abs(pretty(c(stat$freq, -stat$freq))))) +
                       labs(x = NULL, y = "Proportion") +
                       theme(aspect.ratio = plot_ratio,
                             legend.title = element_text(size = text_size),
                             legend.text = element_text(size = text_size),
                             axis.title = element_text(size = text_size),
                             axis.text = element_text(size = text_size),
                       ),
               
               prop_multi = stat %>%
                       mutate(freq = round(freq, 3)) %>%
                       ggplot() + 
                       geom_bar(aes(x = group, y = freq, fill = group), stat = "identity") +
                       geom_text(aes(x = group, y = freq, label = scales::percent(freq)), vjust = -0.5, size = text_size * 0.35) +
                       scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), labels = scales::percent_format()) +
                       facet_wrap(~cluster, ncol = 4, scales = "free") +
                       scale_fill_manual(values = group_colors, name = "Group") +
                       labs(x = NULL, y = "Proportion") + 
                       theme(strip.text.x = element_text(size = text_size)) + thm,
               
               stop("Unknown plot type")
        )
        
        
}

plot_GSEA <- function(gsea_res, pattern = "HALLMARK_", p_cutoff = 0.05, levels) {
        
        gsea_res %>%
                mutate(pathway = str_remove(string = pathway, pattern = pattern)) %>%
                mutate(color = -log10(padj) * sign(NES)) %>%
                mutate(sig = as.factor(ifelse(padj < p_cutoff, 1, 0))) %>%
                ggplot(aes(x = factor(pathway), y = factor(cluster, levels = levels))) + 
                geom_point(aes(size = abs(NES), color = color)) +
                scale_size(name = "Normalized\nEnrichment\nScore Size") +
                scale_color_gradient2(name = bquote(-log[10]~"Adj. p-value"), low = 'dodgerblue1', mid = 'grey', high = 'red', midpoint = 0) +
                coord_flip() +
                geom_point(aes(shape = sig), size = 5.5, stroke = 1) +
                scale_shape_manual(name = "Adj. p-value", values=c(NA, 0), labels = c(paste0("\u2265 ",p_cutoff), paste0("< ",p_cutoff))) +
                theme(axis.title.x = element_blank(),
                      axis.title.y = element_blank())
        
}

plot_measure <- function(dataset, measure, plot_type, group_levels, cluster_levels) {
        
        df <- tibble(group = as.character(dataset$group),
                     cluster = as.character(Idents(dataset)),
                     measure = as.numeric(dataset@meta.data[[measure]]))
        
        thm <- theme(axis.title.x = element_blank(),
                     axis.title.y = element_blank())
        
        switch(plot_type,
               group = ggplot(df, aes(x = factor(group, levels = group_levels), 
                                      y = measure,
                                      fill = factor(group, levels = group_levels))) + 
                       geom_boxplot() +
                       scale_fill_manual(values = get_colors(1:length(group_levels)),
                                         name = "Group") + thm,
               
               cluster = ggplot(df, aes(x = factor(cluster, levels = cluster_levels), 
                                        y = measure,
                                        fill = factor(cluster, levels = cluster_levels))) + 
                       geom_boxplot() +
                       scale_fill_manual(values = get_colors(1:length(cluster_levels)), 
                                         name = "Cluster") + thm,
               
               cluster_group = ggplot(df, aes(x = factor(cluster, levels = cluster_levels), 
                                              y = measure,
                                              fill = factor(group, levels = group_levels))) + 
                       geom_boxplot() +
                       scale_fill_manual(values = get_colors(1:length(group_levels)), 
                                         name = "Group") + thm,
               
               stop("Unknown plot type")
        )
        
}