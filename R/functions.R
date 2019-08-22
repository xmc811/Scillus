
# Load packages and refs for analysis

library(tidyverse)

library(Seurat)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(colormap)
library(ggrepel)
library(reshape2)
library(formattable)
library(magrittr)
library(sctransform)
library(ggpubr)
library(scales)

library(clusterProfiler)
library(org.Mm.eg.db)
library(fgsea)
library(biomaRt)

library(Matrix)
library(fields)
library(KernSmooth)
library(modes)
library(ROCR)
library(DoubletFinder)

library(foreach)
library(doParallel)
library(svglite)

mm_hs <- read_tsv("./refs/mm_hs.txt", col_names = T)
pathways.hallmark <- gmtPathways("./refs/h.all.v6.2.symbols.gmt")

# Helper functions

vlookup <- function(list, table, col_1, col_2) {
        new_list <- table[[col_2]][match(list, table[[col_1]])]
        return(new_list)
}

load_10X_mito_cc <- function(dir = NULL, raw_data = NULL, gcol = 2, proj_name, min_cells = 5, org = "human", cc = T) {
        
        if (!is.null(dir))  raw_data <- Read10X(data.dir = dir, gene.column = gcol)
        
        data <- CreateSeuratObject(counts = raw_data, project = proj_name, min.cells = min_cells)
        
        mito_pattern <- ifelse(org == 'human', "^MT-", "^mt-")
        
        if (org == "mouse") {
                
                s.genes <- vlookup(cc.genes$s.genes, mm_hs, 2, 1)
                g2m.genes <- vlookup(cc.genes$g2m.genes, mm_hs, 2, 1)
                
                s.genes <- s.genes[!is.na(s.genes)]
                g2m.genes <- g2m.genes[!is.na(g2m.genes)]
                
        } else {
                
                s.genes <- cc.genes$s.genes
                g2m.genes <- cc.genes$g2m.genes
        }
        
        data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = mito_pattern)
        
        if (cc) data <- CellCycleScoring(data, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
        
        data$group <- proj_name
        
        return(data)
}

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

filter_norm_10X <- function(dataset, nfeature = 500, mito = 10, nfeatures = 2000) {
        
        expr1 <- FetchData(dataset, vars = "nFeature_RNA")
        expr2 <- FetchData(dataset, vars = "percent.mt")
        
        dataset <- dataset[, which(x = expr1 > nfeature & expr2 < mito)]
        dataset %<>% 
                NormalizeData(verbose = FALSE) %>%
                FindVariableFeatures(selection.method = "vst", nfeatures = nfeatures)
        
}

filter_sctrans_10X <- function(dataset, nfeature = 500, mito = 10, 
                               vars = c("percent.mt","nCount_RNA","S.Score","G2M.Score")) {
        
        expr1 <- FetchData(dataset, vars = "nFeature_RNA")
        expr2 <- FetchData(dataset, vars = "percent.mt")
        
        dataset <- dataset[, which(x = expr1 > nfeature & expr2 < mito)]
        dataset %<>% 
                SCTransform(vars.to.regress = vars, verbose = T)
        
}

find_doublets <- function(dataset, dims = 1:20, ratio = 0.05, resolution = 0.4, txt) {
        
        dataset %<>%
                ScaleData(vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"), verbose = T) %>%
                RunPCA(features = VariableFeatures(dataset)) %>%
                RunUMAP(dims = dims) %>%
                FindNeighbors(dims = dims) %>%
                FindClusters(resolution = resolution)
        
        ## pK Identification
        sweep.res <- paramSweep_v3(dataset, PCs = dims)
        sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
        bcmvn <- find.pK(sweep.stats)
        
        ## Homotypic Doublet Proportion Estimate
        homotypic.prop <- modelHomotypic(dataset$seurat_clusters)
        nExp_poi <- round(ratio*length(Idents(dataset)))
        nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
        
        dataset <- doubletFinder_v3(dataset, PCs = dims, pN = 0.25, pK = 0.1, nExp = nExp_poi.adj, reuse.pANN = F)
        
        barcodes <- names(dataset@active.ident)[dataset[[paste("DF.classifications_0.25_0.1_", as.character(nExp_poi.adj), sep = "")]] == "Doublet"]
        
        paste(barcodes, txt, sep = "")
        
}

analyze_merged <- function(dataset, group.levels,
                           verbose = T, npcs = 50,
                           reduction = "umap", dims = 1:20, nnei = 30, min.dist = 0.3, spread = 1, n.epochs = 500, 
                           k.param = 20,
                           resolution = 0.8) {
        
        dataset$group <- factor(dataset$group, levels = group.levels)
        dataset %<>% 
                RunPCA(npcs = npcs, verbose = verbose) %>%
                RunUMAP(reduction = "pca", dims = dims, n.neighbors = nnei, min.dist = min.dist, spread = spread, n.epochs = n.epochs) %>%
                FindNeighbors(reduction = reduction, dims = 1:2, k.param = k.param) %>%
                FindClusters(resolution = resolution, algorithm = 3)
        
}

get_palette <- function(ncolor, palette = c("Paired", "Set2", "Set1")) {
        
        num <- c()
        for (i in seq(length(palette))) {
                num[i] <- brewer.pal.info[palette[i],][[1]]
        }
        
        ful_pal <- do.call(c, map2(.x = num, .y = palette, .f = brewer.pal))
        
        pal <- ful_pal[1:ncolor]
        
        return(pal)
        
}

get_colors <- function(v, pal = "Paired") {
        
        ncol <- brewer.pal.info[pal,][[1]]
        
        if (sum(!(v %in% 1:ncol)) > 0) {
                stop("Please input a valid numeric vector")
        }
        return(brewer.pal(n = ncol, name = pal)[v])
        
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

# This function is deprecated
find_markers <- function(dataset, use.par = F, ncores = 4, min.diff.pct = -Inf) {
        
        DefaultAssay(object = dataset) <- "RNA"
        
        if (use.par == F) {
                
                markers <- list()
                
                for (i in 0:(length(levels(Idents(dataset)))-1)) {
                        
                        m <- FindConservedMarkers(object = dataset, ident.1 = i, grouping.var = "group", verbose = T, min.diff.pct = min.diff.pct)
                        m %<>%
                                add_column(feature = rownames(m), .before = 1) %>%
                                add_column(cluster = i, .after = 1)
                        
                        markers[[i+1]] <- m
                }
                
        } else {
                
                registerDoParallel(cores = ncores)
                
                markers <- foreach(i = 0:(length(levels(Idents(dataset)))-1)) %dopar% {
                        
                        m <- FindConservedMarkers(object = dataset, ident.1 = i, grouping.var = "group", verbose = T, min.diff.pct = min.diff.pct)
                        m %<>%
                                add_column(feature = rownames(m), .before = 1) %>%
                                add_column(cluster = i, .after = 1)
                        m
                }
        }
        markers <- do.call("rbind", markers)
        
        rownames(markers) <- NULL
        markers <- as_tibble(markers)
        
        return(markers)
}

get_top_genes <- function(dataset, markers, n) {
        
        int_features <- rownames(dataset@assays$integrated@scale.data)
        
        df <- markers %>%
                filter(gene %in% int_features) %>%
                arrange(desc(avg_logFC)) %>%
                group_by(cluster) %>%
                filter(row_number() <= n) %>%
                arrange(cluster)
        
        return(df$gene)
}

plot_heatmap <- function(dataset, markers, nfeatures,
                         cluster_pal = c("Paired", "Set2", "Set1"),
                         group_colors = c('#92c5de','#d6604d')
                         ) {
        
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

rename_cluster <- function(dataset, labels) {
        
        if (length(labels) != length(levels(Idents(dataset)))) {
                stop("Length of new names must be the same with old names.")
        } else {
                current.name <- levels(Idents(dataset))
                new.name <- labels
                
                Idents(dataset) <- plyr::mapvalues(x = Idents(dataset), from = current.name, to = new.name)
                return(dataset)
        }
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

find_diff_genes <- function(dataset, clusters, groups, logfc = 0.25) {
        
        dataset$celltype.group <- paste(Idents(object = dataset), dataset$group, sep = "_")
        Idents(object = dataset) <- "celltype.group"
        
        
        de <- list()
        
        for (i in seq(length(clusters))) {
                
                d <- FindMarkers(dataset, 
                                  ident.1 = paste(clusters[i], groups[2], sep = "_"),
                                  ident.2 = paste(clusters[i], groups[1], sep = "_"),
                                  logfc.threshold = logfc,
                                  assay = "RNA")
                d %<>%
                        add_column(feature = rownames(d), .before = 1) %>%
                        add_column(cluster = clusters[i], .after = 1)
                
                de[[i]] <- as_tibble(d)
        }
        
        de <- do.call("rbind", de)
        
        de

}

test_GSEA <- function(diff, clusters, pathway) {
        
        gsea_res <- list()
        
        for (i in seq(length(clusters))) {
                
                data <- diff %>%
                        filter(cluster == clusters[i]) %>%
                        right_join(mm_hs, by = c("feature" = "mouse")) %>%
                        replace_na(list(avg_logFC = 0)) %>%
                        distinct(human, .keep_all = T) %>%
                        arrange(desc(avg_logFC))
                
                l <- data$avg_logFC
                names(l) <- data$human
                
                res <- fgsea(pathways = pathway, l, minSize = 15, maxSize = 500, nperm = 100000)
                
                res %<>%
                        add_column(cluster = clusters[i], .before = 1)
                
                gsea_res[[i]] <- res
        }
        
        gsea_res <- do.call("rbind", gsea_res)
        
        return(as_tibble(gsea_res))
        
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


add_program_score <- function(dataset, features, org = "human", nbin = 20, ctrl = 10, name){
  
        if(org == "mouse"){
                prog_genes <- vlookup(features, mm_hs, 2, 1)
                prog_genes <- list(prog_genes[!is.na(prog_genes)])
        } else {
                ex_genes <- list(features)
        }
        
        n_genes <- nrow(dataset@assays$integrated@scale.data)
        
        genes_per_bin <- round(n_genes/nbin)
  
        ctrl <- ifelse(ctrl > genes_per_bin, round(genes_per_bin/3), ctrl)
  
        dataset <- AddModuleScore(dataset,
                                  features = prog_genes,
                                  ctrl = ctrl,
                                  nbin = nbin,
                                  name = name)
  
        return(dataset)
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

seurat_to_monocle <- function(dataset, subset = F, clusters = NULL) {
        
        dataset$cluster <- Idents(dataset)
        
        data <- GetAssayData(dataset, assay = 'RNA', slot = 'counts')
        pd <- new('AnnotatedDataFrame', data = dataset@meta.data)
        fd <- new('AnnotatedDataFrame', 
                      data = data.frame(gene_short_name = row.names(data), row.names = rownames(data)))
        
        cds <- newCellDataSet(data,
                                  phenoData = pd,
                                  featureData = fd,
                                  lowerDetectionLimit = 0.5,
                                  expressionFamily = negbinomial.size())
        cds %<>% 
                estimateSizeFactors() %>%
                estimateDispersions() %>%
                detectGenes(min_expr = 0.1)
        
        if (subset) {
                cds <- cds[,row.names(subset(pData(cds), cluster %in% clusters))]
        }
        
        return(cds)

}

analyze_monocle <- function(cds, rev = F) {
        
        expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
        
        diff_test_res <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr = "~group")
        
        ordering_genes <- row.names(subset(diff_test_res, qval < 10E-70))
        
        cds %<>% 
                setOrderingFilter(ordering_genes) %>%
                reduceDimension(max_components = 2, method = 'DDRTree') %>%
                orderCells(reverse = rev)
        
        return(cds)

}


# Test

        
epi_cds <- seurat_to_monocle(dataset = gfp_combined, subset = T, clusters = c('Basal 1','Luminal 1','Basal 2','Luminal 2'))

