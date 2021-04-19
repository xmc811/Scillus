
# These functions are used for visualization


#' Plot the QC metrics of an individual biological sample
#' 
#' @param data_list A list of Seurat objects.
#' @param metrics A string - the name of the QC metrics.
#' @param group_by A string - the grouping variable in the metadata.
#' @param plot_type A string - the type of the plot. Can take values in \code{c("box", "violin", "combined", "density")}.
#' @param pal_setup A dataframe with 2 columns - the \code{RColorBrewer} palette setup to be used; 
#' 
#' Or a character vector of \code{RColorBrewer} palettes. Multiple palettes can be specified, in case of many levels to be colored;
#'  
#' Or a character vector of colors.
#' 
#' @return A plot.
#' @importFrom tibble tibble
#' @importFrom ggplot2 ggplot geom_boxplot geom_violin geom_density aes labs theme_bw waiver
#' @importFrom scales pretty_breaks
#' @importFrom rlang .data
#' @export
#' 

plot_qc <- function(data_list,
                    metrics,
                    group_by = "sample",
                    plot_type = "combined",
                    pal_setup = 'Set2') {
        
        qc <- list()
        
        for (i in seq(length(data_list))) {
                qc[[i]] <- tibble(value = data_list[[i]][[metrics]][[1]],
                                  group = data_list[[i]][[group_by]][[1]])
        }
        
        qc <- do.call(rbind, qc)
        
        if (is.data.frame(pal_setup)) {
            pal <- pal_setup[pal_setup[[1]] == group_by,][[2]]
        } else {
            pal <- pal_setup
        }
        
        colors <- set_colors(pal, length(unique(qc$group)))
        
        p <- ggplot(qc, 
                    mapping = aes(x = .data$group, 
                                  y = .data$value, 
                                  fill = .data$group)) + 
                scale_fill_manual(values = colors) +
                scale_y_continuous(labels = if (metrics == "percent.mt") function(x) paste0(x, "%") else waiver(),
                                   breaks = scales::pretty_breaks(n = 5),
                                   limits = if (!metrics %in% c("S.Score", "G2M.Score")) c(0, NA)) +
                labs(y = metrics) +
                theme_bw() +
                theme(legend.position = "none",
                      axis.title.x = element_blank())
        
        p_den <- ggplot(data = qc,
                        mapping = aes(x = .data$value, 
                                      fill = .data$group)) +
                geom_density(alpha = 0.3, position = "identity") +
                scale_fill_manual(values = colors) +
                labs(x = metrics, y = "Density") +
                theme_bw()
        
        switch(plot_type,
               box = p + geom_boxplot(),
               violin = p + geom_violin(),
               combined = p + geom_violin() + geom_boxplot(alpha = 0.2),
               density = p_den,
               stop("Unknown plot type")
        )
}


#' Single cell visualization - merged and colored by biological sample
#' 
#' @param dataset A Seurat object.
#' @param color_by A string - by which metadata variable the colors will be applied.
#' @param split_by A string - the splitting variable in the metadata.
#' @param pal_setup A dataframe with 2 columns - the \code{RColorBrewer} palette setup to be used;
#' 
#' Or a character vector of \code{RColorBrewer} palettes. Multiple palettes can be specified, in case of many levels to be colored; 
#' 
#' Or a character vector of colors.
#' 
#' @return A plot.
#' @importFrom ggplot2 theme element_rect element_blank facet_wrap scale_color_manual
#' @importFrom stringr str_to_title
#' @importFrom stats as.formula
#' @export
#' 

plot_scdata <- function(dataset, 
                        color_by = "seurat_clusters",
                        split_by = NA,
                        pal_setup = 'Set2') {
        
        if (is.data.frame(pal_setup)) {
                pal <- pal_setup[pal_setup[[1]] == color_by,][[2]]
        } else {
                pal <- pal_setup
        }
    
        split_by <- ifelse(split_by == "No Split", NA, split_by)
    
        color_title <- ifelse(color_by == "seurat_clusters", 
                              "Cluster", 
                              str_to_title(color_by))
        
        df <- dataset@reductions$umap@cell.embeddings %>%
                as.data.frame() %>%
                cbind(dataset@meta.data) %>%
                rownames_to_column(var = "barcode")
        
        colors <- set_colors(pal, length(unique(df[[color_by]])))

        ggplot(df) +
                geom_point(aes(x = .data$UMAP_1,
                               y = .data$UMAP_2,
                               color = .data[[color_by]])) +
                scale_color_manual(values = colors) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      panel.border = element_rect(colour = "black", 
                                                  fill = NA, 
                                                  size = 1, 
                                                  linetype = 1),
                      axis.line = element_blank(),
                      aspect.ratio = 1) +
                labs(x = "UMAP_1", y = "UMAP_2", color = color_title) +
                if (!is.na(split_by)) {facet_wrap(as.formula(paste("~", split_by)))}
}



#' Single cell visualization - gene expression levels
#' 
#' @param dataset A Seurat object.
#' @param measures A character vector - the names of genes or measures to be plotted.
#' @param split_by A string - the splitting variable in the metadata.
#' @param point.size A numeric value - the size of points.
#' 
#' @return A plot.
#' @importFrom ggplot2 element_text xlim ylim ggtitle scale_color_viridis_c
#' @importFrom stats quantile
#' @importFrom purrr map_dbl
#' @export
#' 

plot_measure_dim <- function(dataset, 
                             measures,
                             split_by = NA,
                             point.size = 1) {
        
        l <- get_measure_data(dataset = dataset,
                              measures = measures,
                              return_df = FALSE)
        
        split_by <- ifelse(split_by == "No Split", NA, split_by)
        
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
                        if (!is.na(split_by)) {
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
#' @param markers A dataframe - the dataframe generated by \code{FindAllMarkers} function in \code{Seurat};
#' 
#' Or a character vector - the names of genes can be directly specified.
#' @param sort_var A character vector - the variables to be sorted in the heatmap, should exist in the \code{metadata} of the \code{dataset}.
#' @param n An integer - number of genes to be plotted for each \code{seurat_cluster}. This parameter will not be used if the \code{markers} are directly specified.
#' @param anno_var A character vector - the variables in the annotation, should exist in the \code{metadata} of the \code{dataset}.
#' @param anno_colors A list - the color specification of annotation bars. The length should be the same as \code{anno_var}, with each list element corresponding to each variable. 
#' 
#' For a numeric variable, the element can be the name of a sequential or divergent palette in \code{RColorBrewer}. Please use \code{brewer.pal.info} for more information. 
#' It can also be a character vector of three color values, corresponding to the \code{min}, \code{median}, and \code{max} values of the numeric variable. 
#' 
#' For a categorical variable, the element can be the names of palettes in \code{RColorBrewer}. Multiple palettes can be specified, in case of many levels to be colored. 
#' It can also be a character vector of colors, corresponding to the levels in the categorical variable. Hence, the number of colors specified should be at least the number of levels (unique values) in the categorical variable.
#' @param hm_limit A numeric vector - three numeric values that dictate the lowest limit, midpoint, and highest limit of the color gradient. These three values correspond to the colors specified in the \code{hm_colors}.  
#' @param hm_colors A character vector of length 3 - three colors corresponding to the lowest limit, midpoint, and highest limit specified in \code{hm_colors}.
#' @param row_font_size An integer - the font size of row names.
#' 
#' @return A plot.
#' @importFrom Seurat GetAssayData DefaultAssay
#' @importFrom tibble rownames_to_column
#' @importFrom rlang syms
#' @importFrom circlize colorRamp2
#' @importFrom ggplot2 unit
#' @importFrom ComplexHeatmap Heatmap draw HeatmapAnnotation
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom grid gpar
#' @export
#' 

plot_heatmap <- function(dataset, 
                         markers,
                         sort_var = c('seurat_clusters'),
                         n = 8, 
                         anno_var, 
                         anno_colors,
                         hm_limit = c(-2, 0, 2), 
                         hm_colors = c("#4575b4","white","#d73027"),
                         row_font_size = 12) {
        
        mat <- GetAssayData(object = dataset, assay = DefaultAssay(dataset), slot = "scale.data")
        
        if (is.data.frame(markers)) {
            genes <- get_top_genes(dataset, markers, n)
        } else if (is.character(markers)) {
            genes <- markers
        } else {
            stop('Incorrect input of markers')
        }
        
        
        mat <- mat[match(genes, rownames(mat)),]
        
        anno <- dataset@meta.data %>%
                rownames_to_column(var = "barcode") %>%
                arrange(!!!syms(sort_var))
        
        mat <- t(mat)
        mat <- mat[match(anno$barcode, rownames(mat)),]
        mat <- t(mat)
    
        annos <- list()
        
        for (i in seq_along(1:length(anno_var))) {
            
                err_msg <- paste('Incorrect specification for annotation colors for', anno_var[i])
                
                value <- anno[[anno_var[i]]]
                
                if (is.numeric(value)) {
                    
                    if (all(anno_colors[[i]] %in% rownames(brewer.pal.info)[brewer.pal.info$category != 'qual'])) {
                        
                        n <- brewer.pal.info[anno_colors[[i]],]['maxcolors'][[1]]
                        pal <- brewer.pal(n = n, name = anno_colors[[i]])
                        
                        col_fun <- colorRamp2(c(min(value), stats::median(value), max(value)), 
                                              c(pal[2], pal[(n+1)/2], pal[n-1]))
                        
                    } else if (length(anno_colors[[i]]) == 3 & all(are_colors(anno_colors[[i]]))) {
                        
                        col_fun <- colorRamp2(c(min(value), stats::median(value), max(value)), 
                                              anno_colors[[i]])
                    } else {
                        stop(err_msg)
                    }
                        
                    ha <- HeatmapAnnotation(a = anno[[anno_var[i]]],
                                                col = list(a = col_fun),
                                                border = TRUE,
                                                annotation_label = anno_var[i])
                } else {
                    
                    l <- levels(factor(anno[[anno_var[i]]]))
                    
                    if (all(anno_colors[[i]] %in% rownames(brewer.pal.info))) {
                        
                        col <- set_colors(anno_colors[[i]], length(l))
                        
                    } else if (length(anno_colors[[i]]) >= length(l) & all(are_colors(anno_colors[[i]]))) {
                        
                        col <- anno_colors[[i]]
                        
                    } else {
                        stop(err_msg)
                    }
                        
                    names(col) <- l
                    col <- list(a = col)
                        
                    ha <- HeatmapAnnotation(a = anno[[anno_var[i]]],
                                            col = col,
                                            border = TRUE,
                                            annotation_label = anno_var[i])
                }
                names(ha) <- anno_var[i]
                
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
                      col = colorRamp2(hm_limit, hm_colors),
                      show_column_names = FALSE,
                      row_names_side = "left",
                      row_names_gp = gpar(fontsize = row_font_size),
                      top_annotation = annos)
        
        draw(ht, 
             heatmap_legend_side = "bottom",
             annotation_legend_side = "right")
        
}


#' Plot the basic statistics after cell clustering and labeling
#' 
#' @param dataset A Seurat object.
#' @param plot_type A string - one of \code{"group_count"}, \code{"prop_fill"}, and \code{"prop_multi"}.
#' @param group_by A string - the grouping variable in the metadata.
#' @param pal_setup A dataframe with 2 columns - the \code{RColorBrewer} palette setup to be used; 
#' 
#' or a character vector of \code{RColorBrewer} palettes. Multiple palettes can be specified, in case of many levels to be colored; 
#' 
#' or a character vector of colors.
#' @param plot_ratio A numeric value - the aspect ratio of the plot.
#' @param text_size A numeric value - the text size of the plot.
#' @param tilt_text A logical value - whether the x axis label text is tilted.
#' 
#' @return A plot.
#' @importFrom ggplot2 geom_col geom_text geom_bar scale_fill_manual scale_y_continuous coord_flip facet_wrap expansion
#' @importFrom dplyr mutate summarise n
#' @importFrom forcats fct_rev
#' @importFrom formattable percent
#' @export
#' 

plot_stat <- function(dataset, 
                      plot_type, 
                      group_by = "sample",
                      pal_setup = 'Set2',
                      plot_ratio = 1,
                      text_size = 10,
                      tilt_text = FALSE) {
        
        if (is.data.frame(pal_setup)) {
                pal <- pal_setup[pal_setup[[1]] == group_by,][[2]]
        } else {
                pal <- pal_setup
        }
        
        stat <- tibble::tibble(group = dataset[[group_by]][[1]], 
                               cluster = dataset[['seurat_clusters']][[1]])
        stat %<>%
                group_by(.data$group, 
                         .data$cluster) %>%
                summarise(n = n()) %>%
                mutate(freq = n / sum(n))
        
        ncolors <- if (plot_type == 'prop_fill') {
            length(unique(dataset[['seurat_clusters']][[1]]))
        } else {
            length(unique(dataset[[group_by]][[1]]))
        }
        
        colors <- set_colors(pal, ncolors)
        
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
                       ggplot(aes(x = .data$group, 
                                  y = .data$`sum(n)`)) +
                       geom_col(aes(fill = .data$group)) +
                       geom_text(aes(label = .data$`sum(n)`), 
                                 vjust = -0.5, 
                                 size = text_size * 0.35) +
                       scale_fill_manual(values = colors) + 
                       labs(x = group_by, y = "Number of Cells") + 
                       thm + thm2 + if (tilt_text) {thm3},
               
               prop_fill = ggplot(stat) + 
                       geom_bar(aes(x = .data$group, 
                                    y = .data$freq, 
                                    fill = .data$cluster), 
                                position = "fill", 
                                stat = "identity") +
                       scale_y_continuous(labels = scales::percent) +
                       scale_fill_manual(values = colors, 
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
                       scale_y_continuous(expand = expansion(mult = c(0, 0.1)), 
                                          labels = scales::percent_format()) +
                       facet_wrap(~ cluster, 
                                  ncol = 4, 
                                  scales = "free") +
                       scale_fill_manual(values = colors, 
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
#' @param markers A dataframe - the dataframe generated by \code{FindAllMarkers} function in \code{Seurat}.
#' @param cluster_name A string - the name of the single cell cluster.
#' @param topn An integer - the number of top gene markers in each cluster for the GO enrichment test.
#' @param org A string - the organism that the genes belong to. Can take values from \code{c("human","mouse")}.
#' @param ... Additional arguments to be passed to the function \code{\link[clusterProfiler]{enrichGO}}.
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
#' @param markers A dataframe - the dataframe generated by \code{FindAllMarkers} function in \code{Seurat}.
#' @param org A string - the organism that the genes belong to. Can take values from \code{c("human","mouse")}.
#' @param ... Additional arguments to be passed to the function \code{\link{plot_cluster_go}}.
#' 
#' @return A plot.
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
#' @param gsea_res A dataframe - the results of GSEA generated by \code{test_GSEA}.
#' @param pattern A string - the prefix string pattern to be trimmed of pathway names.
#' @param p_cutoff A numeric value - the p-value cutoff for significance.  
#' @param colors A character vector - three color values, corresponding to the \code{low}, \code{mid}, and \code{high} of the color gradient.
#' 
#' @return A plot.
#' @importFrom stringr str_remove
#' @importFrom ggplot2 geom_point scale_color_gradient2 scale_size scale_shape_manual
#' @export
#' 

plot_GSEA <- function(gsea_res, 
                      pattern = "HALLMARK_", 
                      p_cutoff = 0.05,
                      colors = c('#0570b0','grey','#d7301f')) {
        
        gsea_res %>%
                filter(.data$padj <= p_cutoff) %>%
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
                                      low = colors[1], 
                                      mid = colors[2], 
                                      high = colors[3], 
                                      midpoint = 0) +
                coord_flip() +
                theme_bw() +
                theme(axis.title.x = element_blank(),
                      axis.title.y = element_blank())
}


#' Box plot/Violin plot of gene expressions or meta measures
#' 
#' @param dataset A Seurat object.
#' @param measures A character vector - names of genes or meta measures to plot.
#' @param group_by A string - the grouping variable in the metadata.
#' @param split_by A string - the splitting variable in the metadata.
#' @param pal_setup A dataframe with 2 columns - the \code{RColorBrewer} palette setup to be used; 
#' 
#' Or a character vector of \code{RColorBrewer} palettes. Multiple palettes can be specified, in case of many levels to be colored;
#'  
#' Or a character vector of colors.
#' @param plot_type A string - type of the plot. Should be values in \code{c("box", "violin", "combined")}.
#' 
#' @return A plot.
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom patchwork wrap_plots
#' @export
#' 

plot_measure <- function(dataset, 
                         measures, 
                         group_by,
                         split_by = NA,
                         pal_setup = 'Set2',
                         plot_type = "combined") {
    
    if (is.data.frame(pal_setup)) {
        pal <- pal_setup[pal_setup[[1]] == group_by,][[2]]
    } else {
        pal <- pal_setup
    }
    
    colors <- set_colors(pal, length(unique(dataset[[group_by]][[1]])))
    
    split_by <- ifelse(split_by == "No Split", NA, split_by)
    
    l <- get_measure_data(dataset = dataset,
                          measures = measures,
                          return_df = FALSE)
    
    df <- l[[1]]
    measures <- l[[2]]
    
    p <- list()
    
    a <- ifelse(plot_type == "box", 1, 0.2)
    
    for (i in seq_along(1:length(measures))) {
        
        p[[i]] <- ggplot(df, aes(x = .data[[group_by]],
                                 y = .data[[measures[i]]],
                                 fill = .data[[group_by]])) +
            scale_fill_manual(values = colors) +
            {if (plot_type != "box") geom_violin()} +
            {if (plot_type != "violin") geom_boxplot(alpha = a)} +
            theme(legend.position = "none") +
            labs(x = group_by, y = NULL) +
            {if (!is.na(split_by)) facet_wrap(as.formula(paste0("~", split_by)), nrow = 1)} +
            ggtitle(measures[i])
    }
    patchwork::wrap_plots(p)
}


