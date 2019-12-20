# Thess functions are deprecated but may be useful
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


prop_diverge = stat %>%
        mutate(cluster = fct_rev(.data$cluster)) %>%
        mutate(freq = ifelse(group == group_levels[1], -.data$freq, .data$freq)) %>%
        mutate(freq = round(.data$freq, 3)) %>%
        ggplot() + 
        geom_bar(aes(x = .data$cluster, 
                     y = .data$freq, 
                     fill = .data$group), 
                 stat = "identity") +
        geom_text(aes(x = .data$cluster, 
                      y = .data$freq + 0.03 * sign(.data$freq), 
                      label = percent(abs(.data$freq), 
                                      digits = 1)), 
                  size = text_size * 0.35) +
        coord_flip() +
        scale_fill_manual(values = group_colors, name = "Group") +
        scale_y_continuous(breaks = pretty(c(stat$freq, -stat$freq)),
                           labels = scales::percent(abs(pretty(c(stat$freq, -stat$freq))))) +
        labs(x = NULL, 
             y = "Proportion") +
        theme(aspect.ratio = plot_ratio,
              legend.title = element_text(size = text_size),
              legend.text = element_text(size = text_size),
              axis.title = element_text(size = text_size),
              axis.text = element_text(size = text_size),
        )
