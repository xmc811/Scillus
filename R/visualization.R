
# These functions are used for visualization


#' Plot the QC metrics of an individual biological sample
#' 
#' @param data_list A list of Seurat objects
#' @param metrics A string - the name of the QC metrics
#' @param plot_type A string - the type of the plot. Default value is \code{"combined"}.
#' @param palette A string - the \code{RColorBrewer} palette to be used. Default value is \code{"Set2"}.
#' 
#' @return A plot.
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot geom_boxplot geom_violin aes labs theme_bw waiver
#' @importFrom scales pretty_breaks
#' @importFrom rlang .data
#' @export
#' 

plot_qc <- function(data_list, 
                    metrics, 
                    plot_type = "combined",
                    palette = "Set2") {
        
        qc <- list()
        
        for (i in seq(length(data_list))) {
                qc[[i]] <- tibble(value = data_list[[i]][[metrics]][[1]],
                                  sample = data_list[[i]]@project.name)
        }
        
        qc <- do.call(rbind, qc)
        
        p <- ggplot(qc, 
                    mapping = aes(x = .data$sample, 
                                  y = .data$value, 
                                  fill = sample)) + 
                scale_fill_manual(values = get_spectrum(n = length(data_list), 
                                                        palette = palette)) +
                scale_y_continuous(labels = if (metrics == "percent.mt") function(x) paste0(x, "%") else waiver(),
                                   breaks = scales::pretty_breaks(n = 5),
                                   limits = if (!metrics %in% c("S.Score", "G2M.Score")) c(0, NA)) +
                labs(y = metrics) +
                theme_bw() +
                theme(legend.position = "none",
                      axis.text = element_text(size = 12),
                      axis.title.y = element_text(size = 12),
                      axis.title.x = element_blank())
        
        switch(plot_type,
               box = p + geom_boxplot(),
               violin = p + geom_violin(),
               combined = p + geom_violin() + geom_boxplot(alpha = 0.2),
               stop("Unknown plot type")
        )
}


#' Single cell visualization - merged and colored by biological sample
#' 
#' @param dataset A Seurat object.
#' @param colors A string vector - the colors used for samples. Default value is \code{NULL}.
#' @param color_by A string - by which metadata the colors will be applied. Default value is \code{"seurat_clusters"}.
#' @param split A string - by which metadata the plot will be split. Default value is \code{NULL}.
#' @param ... Arguments passed to \code{DimPlot}.
#' 
#' @return A plot.
#' @importFrom Seurat DimPlot
#' @importFrom ggplot2 scale_color_manual theme element_rect element_blank
#' @export
#' 

plot_scdata <- function(dataset, 
                        color_by = "seurat_clusters", 
                        colors = NULL, 
                        split = NULL, 
                        ...) {
        
        dataset$sample <- factor(dataset$sample)
        
        if (color_by == "sample") {
                levels <- levels(dataset$sample)
        } else {
                levels <- levels(dataset$seurat_clusters)
        }
        
        legend.title <- ifelse(color_by == "sample", "Sample", "Cluster")
        
        if (is.null(colors)) {
                if (color_by == "sample") {
                        colors <- get_spectrum(length(levels))
                } else {
                        colors <- get_palette(length(levels))
                }
                
        }
        
        thm <- theme(panel.border = element_rect(colour = "black", 
                                                 fill = NA, 
                                                 size = 1, 
                                                 linetype = 1),
                     axis.line=element_blank(),
                     aspect.ratio = 1)
        
        if (is.null(split)) {
                p <- DimPlot(object = dataset, 
                             reduction = "umap", 
                             group.by = color_by, ...)
                
                p + scale_color_manual(values = colors, 
                                       name = legend.title, 
                                       labels = levels) + thm
                
        } else {
                p <- DimPlot(object = dataset, 
                             reduction = "umap", 
                             split.by = split, 
                             group.by = color_by, 
                             ncol = ceiling(sqrt(length(levels(dataset$sample)))), ...)
                
                p  + scale_color_manual(values = colors, 
                                        name = legend.title, 
                                        labels = levels) + 
                        thm + 
                        theme(strip.text.x = element_text(face = "plain", 
                                                          vjust = 1))
        }
        
}



#' Single cell visualization - gene expression levels
#' 
#' @param dataset A Seurat object.
#' @param features A string vector -
#' @param plot_type A string -
#' @param meta A logical value -
#' 
#' @return A plot.
#' @importFrom Seurat FeaturePlot
#' @importFrom ggplot2 element_text xlim ylim ggtitle scale_color_viridis_c
#' @importFrom grid grobTree rectGrob textGrob gpar viewport
#' @importFrom stats quantile
#' @importFrom colormap colormap
#' @importFrom purrr map_dbl
#' @importFrom tidyr nest
#' @export
#' 

plot_measure_cluster <- function(dataset, 
                                 features, 
                                 plot_type = "none", 
                                 meta = FALSE) {
        
        dataset$group <- factor(dataset$group)
        Idents(dataset) <- factor(Idents(dataset))
        
        df <- if (meta) get_meta_data(dataset, features) else get_gene_data(dataset, features)
        
        df$group <- factor(df$group, levels = levels(dataset$group))
        df$cluster <- factor(df$cluster, levels = levels(Idents(dataset)))
        
        lst <- df %>% group_by(feature) %>% nest()
        
        plt_func <- function(df, title, type) {
                ggplot(df, aes(x = x_cor, y = y_cor, color = value)) + 
                        geom_point(size = 0.2) +
                        scale_color_viridis_c(option = "A",
                                              name = "",
                                              direction = -1, 
                                              limits = c(quantile(df$value, probs = 0.1), max(df$value)), 
                                              oob = scales::squish) +
                        ggtitle(title) +
                        theme(panel.background = element_rect(fill = 'white', colour = 'black'),
                              panel.border = element_rect(colour = "black", fill = NA, size = 1, linetype = 1),
                              aspect.ratio = 1,
                              plot.title = element_text(hjust = 0.5),
                              axis.title = element_blank()
                        ) +
                        if (type == "group") {
                                facet_wrap( ~ group)
                        } else if (type == "cluster") {
                                facet_wrap( ~ cluster)
                        } else {}
                
        }
        
        l <- map2(.x = lst$data, .y = lst$feature, type = plot_type, .f = plt_func)
        
        grid.arrange(grobs = l, ncol = ceiling(sqrt(length(l))))
        
}


#' Plot the heatmap of single cell dataset
#' 
#' @param dataset A Seurat object.
#' @param markers A tibble -
#' @param sort_var A string vector -
#' @param n An integer -
#' @param anno_var A string vector -
#' @param anno_colors A string vector -
#' 
#' @return A plot.
#' @importFrom Seurat GetAssayData
#' @importFrom tibble rownames_to_column
#' @importFrom rlang syms
#' @importFrom circlize colorRamp2
#' @importFrom ggplot2 unit
#' @importFrom ComplexHeatmap Heatmap draw HeatmapAnnotation
#' @export
#' 

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
                        
                        col_fun <- colorRamp2(c(min(value), stats::median(value), max(value)), 
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


#' Plot the basic statistics after cell clustering and labeling
#' 
#' @param dataset A Seurat object.
#' @param plot_type A string -
#' @param cluster_colors A string vector -
#' @param sample_colors A string vector -
#' @param plot_ratio A double -
#' @param text_size A double -
#' @param tilt_text A logical -
#' 
#' @return A plot.
#' @importFrom ggplot2 geom_col geom_text geom_bar scale_fill_manual scale_y_continuous coord_flip facet_wrap expand_scale
#' @importFrom dplyr mutate summarise n
#' @importFrom forcats fct_rev
#' @importFrom formattable percent
#' @export
#' 

plot_stat <- function(dataset, 
                      plot_type, 
                      cluster_colors = NULL,
                      sample_colors = NULL,
                      plot_ratio = 1,
                      text_size = 10,
                      tilt_text = FALSE) {
        
        dataset$sample <- factor(dataset$sample)
        
        sample_levels <- levels(dataset$sample)
        cluster_levels <- levels(dataset$seurat_clusters)
        
        if (is.null(sample_colors)) sample_colors <- get_spectrum(length(sample_levels))
        if (is.null(cluster_colors)) cluster_colors <- get_palette(length(cluster_levels))
        
        stat <- as_tibble(cbind(sample = as.character(dataset$sample), cluster = as.character(Idents(dataset))))
        stat %<>%
                mutate(sample = factor(.data$sample, levels = sample_levels),
                       cluster = factor(.data$cluster, levels = cluster_levels)) %>%
                group_by(.data$sample, .data$cluster) %>%
                summarise(n = n()) %>%
                mutate(freq = n / sum(n))
        
        
        thm <- theme(aspect.ratio = plot_ratio,
                     legend.title = element_text(size = text_size),
                     legend.text = element_text(size = text_size),
                     axis.title = element_text(size = text_size),
                     axis.text = element_text(size = text_size),
                     axis.title.x = element_blank()
        ) + theme_bw()
        thm2 <- theme(legend.position = "none")
        thm3 <- theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
        
        switch(plot_type,
               sample_count = stat %>%
                       group_by(.data$sample) %>%
                       summarise(sum(n)) %>%
                       ggplot() +
                       geom_col(aes(x = .data$sample, 
                                    y = .data$`sum(n)`, 
                                    fill = .data$sample)) +
                       geom_text(aes(x = .data$sample, 
                                     y = .data$`sum(n)`, 
                                     label = .data$`sum(n)`), 
                                 vjust = -0.5, 
                                 size = text_size * 0.35) +
                       scale_fill_manual(values = sample_colors) + 
                       labs(y = "Number of Cells") + 
                       thm + thm2 + if (tilt_text) {thm3},
               
               cluster_count = stat %>%
                       group_by(.data$cluster) %>%
                       summarise(sum(n)) %>%
                       ggplot() +
                       geom_col(aes(x = .data$cluster, 
                                    y = .data$`sum(n)`, 
                                    fill = .data$cluster)) +
                       geom_text(aes(x = .data$cluster, 
                                     y = .data$`sum(n)`, 
                                     label = .data$`sum(n)`), 
                                 vjust = -0.5, 
                                 size = text_size * 0.35) +
                       scale_fill_manual(values = cluster_colors, name = "Cluster") + 
                       labs(y = "Number of Cells") + 
                       thm + thm2 + if (tilt_text) {thm3},
               
               prop_fill = ggplot(stat) + 
                       geom_bar(aes(x = .data$sample, 
                                    y = .data$freq, 
                                    fill = .data$cluster), 
                                position = "fill", 
                                stat = "identity") +
                       scale_y_continuous(labels = scales::percent) +
                       scale_fill_manual(values = cluster_colors, name = "Cluster") +
                       labs(y = "Proportion") + 
                       thm + if (tilt_text) {thm3},
               
               prop_diverge = stat %>%
                       mutate(cluster = fct_rev(.data$cluster)) %>%
                       mutate(freq = ifelse(sample == sample_levels[1], -.data$freq, .data$freq)) %>%
                       mutate(freq = round(.data$freq, 3)) %>%
                       ggplot() + 
                       geom_bar(aes(x = .data$cluster, 
                                    y = .data$freq, 
                                    fill = .data$sample), 
                                stat = "identity") +
                       geom_text(aes(x = .data$cluster, 
                                     y = .data$freq + 0.03 * sign(.data$freq), 
                                     label = percent(abs(.data$freq), digits = 1)), 
                                 size = text_size * 0.35) +
                       coord_flip() +
                       scale_fill_manual(values = sample_colors, name = "Sample") +
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
                       mutate(freq = round(.data$freq, 3)) %>%
                       ggplot() + 
                       geom_bar(aes(x = .data$sample,
                                    y = .data$freq, 
                                    fill = .data$sample), 
                                stat = "identity") +
                       geom_text(aes(x = .data$sample, 
                                     y = .data$freq, 
                                     label = scales::percent(.data$freq)), 
                                 vjust = -0.5, 
                                 size = text_size * 0.35) +
                       scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), labels = scales::percent_format()) +
                       facet_wrap(~cluster, ncol = 4, scales = "free") +
                       scale_fill_manual(values = sample_colors, name = "Sample") +
                       labs(x = NULL, y = "Proportion") + 
                       theme(strip.text.x = element_text(size = text_size)) + 
                       thm + thm2 + if (tilt_text) {thm3},
               
               stop("Unknown plot type")
        )
}


#' plot the GO enrichment analysis of a cluster
#' 
#' @param markers A tibble -
#' @param cluster_name A string -
#' @param topn An integer -
#' @param org A string -
#' @param ... Additional arguments to be passed to the function \code{\link{enrichGO}}.
#' 
#' @return A plot.
#' @importFrom clusterProfiler enrichGO
#' @importFrom utils head
#' @importFrom dplyr pull
#' @importFrom forcats fct_reorder
#' @importFrom ggplot2 facet_grid geom_hline xlab ylab scale_fill_gradient2
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @export
#' 

plot_cluster_go <- function (markers, cluster_name, topn = 100, org, ...) {
        
        gene_list <- markers %>% 
                filter(.data$cluster == cluster_name) %>% 
                arrange(.data$p_val_adj) %>% 
                head(topn) %>% 
                pull(.data$gene)
        
        db <- if (org == "human") org.Hs.eg.db::org.Hs.eg.db else org.Mm.eg.db::org.Mm.eg.db
        
        res <- enrichGO(gene = gene_list, 
                        OrgDb = db, 
                        keyType = "SYMBOL", ...)
        
        df <- as_tibble(res@result) %>% 
                arrange(.data$p.adjust) %>% 
                head(10) %>% 
                mutate(cluster = cluster_name) %>% 
                mutate(Description = stringr::str_to_title(.data$Description)) %>% 
                mutate(Description = fct_reorder(.data$Description, dplyr::desc(.data$p.adjust)))
        
        ggplot(df, 
               mapping = aes(x = .data$Description, 
                             y = -log10(.data$p.adjust))) + 
                geom_bar(aes(fill = .data$Count), 
                         stat = "identity") + 
                scale_fill_gradient2("Gene Count", 
                                     low = "lightgrey", 
                                     mid = "#feb24c", 
                                     high = "#bd0026") + 
                coord_flip() + 
                geom_hline(yintercept = -log10(0.05), 
                           linetype = "dashed") + 
                xlab("Gene Ontology") + 
                ylab(bquote("-log"[10] ~ " adjusted p-value")) + facet_grid(. ~ .data$cluster)
}


#' plot the GO enrichment analysis of all clusters of a dataset
#' 
#' @param markers A tibble -
#' @param org A string -
#' @param ... Additional arguments to be passed to the function \code{\link{plot_cluster_go}}.
#' 
#' @return A plot.
#' @importFrom gridExtra grid.arrange
#' @export
#' 

plot_all_cluster_go <- function (markers, org = "human", ...) {
        
        lst <- list()
        clusters <- levels(markers$cluster)
        lst <- purrr::map(.x = clusters, .f = plot_cluster_go, markers = markers, org = org, ...)
        do.call("grid.arrange", c(lst, ncol = floor(sqrt(length(lst)))))
        
}


#' plot the GSVA scores of given gene signatures
#' 
#' @param dataset A Seurat object
#' @param gene_list A list -
#' @param pattern A string -
#' @param ... Additional arguments to be passed to the function \code{\link{gsva}}.
#' 
#' @return A plot.
#' @importFrom Seurat AverageExpression
#' @importFrom stats dist hclust
#' @importFrom GSVA gsva
#' @importFrom ggplot2 geom_tile
#' @export
#' 

plot_GSVA <- function(dataset, gene_list, pattern = "HALLMARK_", ...) {
        
        names(gene_list) <- str_remove(names(gene_list), pattern = pattern)
        
        mtx <- as.matrix(AverageExpression(dataset, assays = "RNA")[[1]])
        
        res <- gsva(expr = mtx, gset.idx.list = gene_list, ...)
        
        d <- dist(res)
        
        lvl <- hclust(d)$labels[hclust(d)$order]
        
        res2 <- res %>% 
                as_tibble(rownames = "pathway") %>%
                melt(id = "pathway")
        
        ggplot(res2, aes(x = factor(pathway, levels = rev(lvl)), y = variable)) + 
                geom_tile(aes(fill = value), color = "black") +
                scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0) +
                coord_flip() +
                labs(x = "Gene Signatures", y = "Clusters", fill = "GSVA score")
        
}


#' plot the results of GSEA
#' 
#' @param gsea_res A tibble -
#' @param pattern A string -
#' @param p_cutoff A double -
#' @param levels A string vector -
#' 
#' @return A plot.
#' @importFrom stringr str_remove
#' @importFrom ggplot2 geom_point scale_color_gradient2 scale_size scale_shape_manual
#' @export
#' 

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


#' Box plot/Violin plot of gene expressions or meta measures
#' 
#' @param dataset A Seurat object.
#' @param features A string vector -
#' @param plot_type A string -
#' @param meta A logical value -
#' @param group_colors A string vector -
#' @param cluster_colors A string vector -
#' @param show A string -
#' 
#' @return A plot.
#' @export
#' 

plot_measure <- function(dataset, 
                         features, 
                         plot_type, 
                         show = "combined",
                         meta = FALSE,
                         group_colors = NULL, 
                         cluster_colors = NULL) {
        
        dataset$group <- factor(dataset$group)
        Idents(dataset) <- factor(Idents(dataset))
        
        group_levels <- levels(dataset$group)
        cluster_levels <- levels(Idents(dataset))
        
        if (is.null(group_colors)) group_colors <- get_spectrum(length(group_levels))
        if (is.null(cluster_colors)) cluster_colors <- get_palette(length(cluster_levels))
        
        df <- if (meta) get_meta_data(dataset, features) else get_gene_data(dataset, features)
        
        df$feature <- factor(df$feature, levels = features)
        
        n <- ceiling(sqrt(length(unique(df$feature))))
        
        thm <- theme(axis.title.y = element_blank())
        thm2 <- theme(legend.position = "none")
        
        a <- ifelse(show == "box", 1, 0.2)
        
        switch(plot_type,
               group = ggplot(df, aes(x = factor(group, levels = group_levels), 
                                      y = value,
                                      fill = factor(group, levels = group_levels))) + 
                       {if (show != "box") geom_violin()} +
                       {if (show != "violin") geom_boxplot(alpha = a)} +
                       xlab("Sample") +
                       scale_fill_manual(values = group_colors) + thm + thm2 +
                       facet_wrap( ~ feature, scales = "free", ncol = n),
               
               cluster = ggplot(df, aes(x = factor(cluster, levels = cluster_levels), 
                                        y = value,
                                        fill = factor(cluster, levels = cluster_levels))) + 
                       {if (show != "box") geom_violin()} + 
                       {if (show != "violin") geom_boxplot(alpha = a)} +
                       xlab("Cluster") +
                       scale_fill_manual(values = cluster_colors) + thm + thm2 +
                       facet_wrap( ~ feature, scales = "free", ncol = n),
               
               cluster_group = ggplot(df, aes(x = factor(cluster, levels = cluster_levels), 
                                              y = value,
                                              fill = factor(group, levels = group_levels))) + 
                       {if (show != "box") geom_violin(position = position_dodge(width = 0.8))} + 
                       {if (show != "violin") geom_boxplot(alpha = a, width = 0.7, position = position_dodge(width = 0.8))} +
                       xlab("Cluster") +
                       scale_fill_manual(values = group_colors, name = "Sample") + thm +
                       facet_wrap( ~ feature, scales = "free", ncol = n),
               
               stop("Unknown plot type")
        )
}
