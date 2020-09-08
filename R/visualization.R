
# These functions are used for visualization


#' Plot the QC metrics of an individual biological sample
#' 
#' @param data_list A list of Seurat objects
#' @param metrics A string - the name of the QC metrics
#' @param group_by A string - the grouping variable for plot in the metadata. Default value is \code{"sample"}.
#' @param plot_type A string - the type of the plot. Default value is \code{"combined"}.
#' @param pal_setup A dataframe with 2 columns - the \code{RColorBrewer} palette setup to be used. Default value is \code{NULL}.
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
                    group_by = "sample",
                    plot_type = "combined",
                    pal_setup = NULL) {
        
        qc <- list()
        
        for (i in seq(length(data_list))) {
                qc[[i]] <- tibble(value = data_list[[i]][[metrics]][[1]],
                                  group = data_list[[i]][[group_by]][[1]])
        }
        
        qc <- do.call(rbind, qc)
        
        if(!is.null(pal_setup)) {
                palette <- pal_setup[pal_setup[[1]] == group_by,][[2]]
        } else {
                palette <- "Set2"
        }
        
        p <- ggplot(qc, 
                    mapping = aes(x = .data$group, 
                                  y = .data$value, 
                                  fill = .data$group)) + 
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
#' @param color_by A string - by which metadata the colors will be applied. Default value is \code{"seurat_clusters"}.
#' @param split_by A string - by which metadata the plot will be split. Default value is \code{NULL}.
#' @param pal_setup A dataframe with 2 columns or a string - the \code{RColorBrewer} palette setup to be used. Default value is \code{NULL}.
#' 
#' @return A plot.
#' @importFrom ggplot2 theme element_rect element_blank facet_wrap scale_color_brewer
#' @importFrom stringr str_to_title
#' @importFrom stats as.formula
#' @export
#' 

plot_scdata <- function(dataset, 
                        color_by = "seurat_clusters",
                        split_by = NULL,
                        pal_setup = NULL) {
        
        if (is.data.frame(pal_setup)) {
                pal <- pal_setup[pal_setup[[1]] == color_by,][[2]]
        } else if (is.character(pal_setup)){
                pal <- pal_setup
        } else {
                pal <- "Set2"
        }
    
        color_title <- ifelse(color_by == "seurat_clusters", 
                              "Cluster", 
                              str_to_title(color_by))
        
        df <- dataset@reductions$umap@cell.embeddings %>%
                as.data.frame() %>%
                cbind(dataset@meta.data) %>%
                rownames_to_column(var = "barcode")

        ggplot(df) +
                geom_point(aes(x = .data$UMAP_1,
                               y = .data$UMAP_2,
                               color = .data[[color_by]])) +
                scale_color_brewer(palette = pal) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.text = element_text(size = 12),
                      axis.title = element_text(size = 12),
                      panel.border = element_rect(colour = "black", 
                                                  fill = NA, 
                                                  size = 1, 
                                                  linetype = 1),
                      axis.line = element_blank(),
                      aspect.ratio = 1) +
                labs(x = "UMAP_1", y = "UMAP_2", color = color_title) +
                if (!is.null(split_by)) {facet_wrap(as.formula(paste("~", split_by)))}
}



#' Single cell visualization - gene expression levels
#' 
#' @param dataset A Seurat object.
#' @param measures A string vector -
#' @param split_by A string -
#' @param point.size A double -
#' 
#' @return A plot.
#' @importFrom ggplot2 element_text xlim ylim ggtitle scale_color_viridis_c
#' @importFrom stats quantile
#' @importFrom purrr map_dbl
#' @export
#' 

plot_measure_dim <- function(dataset, 
                             measures,
                             split_by = NULL,
                             point.size = 1) {
        
        l <- get_measure_data(dataset = dataset,
                              measures = measures,
                              return_df = FALSE)
        
        df <- l[[1]]
        measures <- l[[2]]
        
        p <- list()
        
        for (i in seq_along(1:length(measures))) {
    
                p[[i]] <- ggplot(df) +
                        geom_point(aes(x = .data$UMAP_1,
                                       y = .data$UMAP_2,
                                       color = .data[[measures[i]]]),
                                   size = point.size) +
                        scale_color_viridis_c(option = "A",
                                              name = "",
                                              direction = -1, 
                                              limits = c(quantile(.data[[measures[i]]], 
                                                                  probs = 0.1), 
                                                         quantile(.data[[measures[i]]], 
                                                                  probs = 0.9)), 
                                              oob = scales::squish) +
                        theme(panel.grid.major = element_blank(),
                              panel.grid.minor = element_blank(),
                              panel.background = element_blank(),
                              axis.text = element_text(size = 12),
                              axis.title = element_text(size = 12),
                              panel.border = element_rect(colour = "black", 
                                                          fill = NA, 
                                                          size = 1, 
                                                          linetype = 1),
                              axis.line = element_blank()) +
                        labs(x = "UMAP_1", y = "UMAP_2") +
                        ggtitle(measures[i]) +
                        if (!is.null(split_by)) {
                                facet_wrap(as.formula(paste("~", split_by)))
                        } else {
                                theme(aspect.ratio = 1)
                        }
                
        }
        patchwork::wrap_plots(p)
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
#' @param plot_type A string - one of \code{"group_count"}, \code{"cluster_count"}, \code{"prop_fill"}, and \code{"prop_multi"}.
#' @param group_by A string - the grouping variable for plot in the metadata. Default value is \code{"sample"}.
#' @param pal_setup A dataframe with 2 columns - the \code{RColorBrewer} palette setup to be used. Default value is \code{NULL}.
#' @param plot_ratio A double - the aspect ratio of the plot. Default value is \code{1}.
#' @param text_size A double - the text size of the plot. Default value is \code{10}.
#' @param tilt_text A logical - whether the x axis label text is tilted. Default value is \code{FALSE}.
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
                      group_by = "sample",
                      pal_setup = NULL,
                      plot_ratio = 1,
                      text_size = 10,
                      tilt_text = FALSE) {
        
        if (!is.null(pal_setup)) {
                group_palette <- pal_setup[pal_setup[[1]] == group_by,][[2]]
                cluster_palette <- pal_setup[pal_setup[[1]] == "seurat_clusters",][[2]]
        } else {
                group_palette <- "Set2"
                cluster_palette <- "Paired"
        }
              
        group_colors <- get_spectrum(n = length(unique(dataset[[group_by]][[1]])), 
                                     palette = group_palette)
        
        cluster_colors <- get_spectrum(n = length(unique(Idents(dataset))), 
                                       palette = cluster_palette)
        
        stat <- tibble::tibble(group = dataset[[group_by]][[1]], 
                               cluster = Idents(dataset))
        stat %<>%
                group_by(.data$group, 
                         .data$cluster) %>%
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
        thm3 <- theme(axis.text.x = element_text(angle = 45, 
                                                 vjust = 0.5))
        
        switch(plot_type,
               group_count = stat %>%
                       group_by(.data$group) %>%
                       summarise(sum(n)) %>%
                       ggplot() +
                       geom_col(aes(x = .data$group, 
                                    y = .data$`sum(n)`, 
                                    fill = .data$group)) +
                       geom_text(aes(x = .data$group, 
                                     y = .data$`sum(n)`, 
                                     label = .data$`sum(n)`), 
                                 vjust = -0.5, 
                                 size = text_size * 0.35) +
                       scale_fill_manual(values = group_colors) + 
                       labs(x = group_by, y = "Number of Cells") + 
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
                       scale_fill_manual(values = cluster_colors, 
                                         name = "Cluster") + 
                       labs(x = "Cluster", y = "Number of Cells") + 
                       thm + thm2 + if (tilt_text) {thm3},
               
               prop_fill = ggplot(stat) + 
                       geom_bar(aes(x = .data$group, 
                                    y = .data$freq, 
                                    fill = .data$cluster), 
                                position = "fill", 
                                stat = "identity") +
                       scale_y_continuous(labels = scales::percent) +
                       scale_fill_manual(values = cluster_colors, 
                                         name = "Cluster") +
                       labs(x = group_by, y = "Proportion") + 
                       thm + if (tilt_text) {thm3},
               
               prop_multi = stat %>%
                       mutate(freq = round(.data$freq, 3)) %>%
                       ggplot() + 
                       geom_bar(aes(x = .data$group,
                                    y = .data$freq, 
                                    fill = .data$group), 
                                stat = "identity") +
                       geom_text(aes(x = .data$group, 
                                     y = .data$freq, 
                                     label = scales::percent(.data$freq)), 
                                 vjust = -0.5, 
                                 size = text_size * 0.35) +
                       scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), 
                                          labels = scales::percent_format()) +
                       facet_wrap(~ cluster, 
                                  ncol = 4, 
                                  scales = "free") +
                       scale_fill_manual(values = group_colors, 
                                         name = "Group") +
                       labs(x = NULL, 
                            y = "Proportion") + 
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
#' @param ... Additional arguments to be passed to the function \code{\link[clusterProfiler]{enrichGO}}
#' 
#' @return A plot.
#' @importFrom utils head
#' @importFrom dplyr pull
#' @importFrom forcats fct_reorder
#' @importFrom ggplot2 facet_grid geom_hline xlab ylab scale_fill_gradient2
#' @export
#' 

plot_cluster_go <- function (markers, cluster_name, topn = 100, org, ...) {
        
        pkg_name <- ifelse(org == "human", "org.Hs.eg.db", "org.Mm.eg.db")
        
        if (!requireNamespace(pkg_name, quietly = TRUE)) {
                stop(paste("Package", pkg_name, "needed for this function to work. Please install it."),
                     call. = FALSE)
        }
        
        if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
                stop(paste("Package \"clusterProfiler\" needed for this function to work. Please install it."),
                     call. = FALSE)
        }
        
        gene_list <- markers %>% 
                filter(.data$cluster == cluster_name) %>% 
                arrange(.data$p_val_adj) %>% 
                head(topn) %>% 
                pull(.data$gene)
        
        db <- if(org == "human") org.Hs.eg.db::org.Hs.eg.db else org.Mm.eg.db::org.Mm.eg.db
        
        res <- clusterProfiler::enrichGO(gene = gene_list, 
                                         OrgDb = db, 
                                         keyType = "SYMBOL", ...)
        
        df <- as_tibble(res@result) %>% 
                arrange(.data$p.adjust) %>% 
                head(10) %>% 
                mutate(cluster = cluster_name) %>% 
                mutate(Description = stringr::str_to_title(.data$Description)) %>% 
                mutate(Description = fct_reorder(.data$Description, 
                                                 dplyr::desc(.data$p.adjust)))
        
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
                ylab(bquote("-log"[10] ~ " adjusted p-value")) +
                theme_bw() +
                theme(axis.text = element_text(size = 10),
                      axis.title = element_text(size = 12)) +
                ggtitle(cluster_name)
}


#' plot the GO enrichment analysis of all clusters of a dataset
#' 
#' @param markers A tibble -
#' @param org A string -
#' @param ... Additional arguments to be passed to the function \code{\link{plot_cluster_go}}.
#' 
#' @return A plot
#' @export
#' 

plot_all_cluster_go <- function (markers, org = "human", ...) {
        
        lst <- list()
        clusters <- levels(markers$cluster)
        lst <- purrr::map(.x = clusters, 
                          .f = plot_cluster_go, 
                          markers = markers, 
                          org = org, ...)
        
        patchwork::wrap_plots(lst, ncol = 2)
}


#' plot the results of GSEA
#' 
#' @param gsea_res A tibble -
#' @param pattern A string -
#' @param p_cutoff A double -
#' 
#' @return A plot.
#' @importFrom stringr str_remove
#' @importFrom ggplot2 geom_point scale_color_gradient2 scale_size scale_shape_manual
#' @export
#' 

plot_GSEA <- function(gsea_res, 
                      pattern = "HALLMARK_", 
                      p_cutoff = 0.05) {
        
        gsea_res %>%
                filter(.data$padj <= 0.05) %>%
                mutate(pathway = str_remove(string = .data$pathway, 
                                            pattern = pattern)) %>%
                mutate(color = -log10(.data$padj) * sign(.data$NES)) %>%
                ggplot() + 
                geom_point(aes(x = factor(.data$pathway), 
                               y = factor(.data$cluster),
                               size = abs(.data$NES), 
                               color = .data$color)) +
                scale_size(name = "Normalized\nEnrichment\nScore Size") +
                scale_color_gradient2(name = bquote(-log[10]~"Adj. p-value"), 
                                      low = '#0570b0', 
                                      mid = 'grey', 
                                      high = '#d7301f', 
                                      midpoint = 0) +
                coord_flip() +
                theme_bw() +
                theme(axis.title.x = element_blank(),
                      axis.title.y = element_blank())
}


#' Box plot/Violin plot of gene expressions or meta measures
#' 
#' @param dataset A Seurat object.
#' @param measures A string vector -
#' @param group_by A string -
#' @param pal_setup A dataframe with 2 columns - the \code{RColorBrewer} palette setup to be used. Default value is \code{NULL}.
#' @param show A string -
#' 
#' @return A plot.
#' @importFrom ggplot2 scale_fill_brewer
#' @importFrom patchwork wrap_plots
#' @export
#' 

plot_measure <- function(dataset, 
                         measures, 
                         group_by,
                         pal_setup = NULL,
                         show = "combined") {
        
        if (!is.null(pal_setup)) {
                pal <- pal_setup[pal_setup[[1]] == group_by,][[2]]
        } else {
                pal <- "Set2"
        }
        
        l <- get_measure_data(dataset = dataset,
                              measures = measures,
                              return_df = FALSE)
        
        df <- l[[1]]
        measures <- l[[2]]
        
        p <- list()
        
        a <- ifelse(show == "box", 1, 0.2)
        
        for (i in seq_along(1:length(measures))) {
                
                p[[i]] <- ggplot(df, aes(x = .data[[group_by]],
                                         y = .data[[measures[i]]],
                                         fill = .data[[group_by]])) +
                        scale_fill_brewer(palette = pal) +
                        {if (show != "box") geom_violin()} +
                        {if (show != "violin") geom_boxplot(alpha = a)} +
                        theme(legend.position = "none") +
                        labs(x = group_by, y = NULL) +
                        ggtitle(measures[i])
        }
        patchwork::wrap_plots(p)
}
